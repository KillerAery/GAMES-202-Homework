#ifdef GL_ES
precision mediump float;
#endif

// Phong related variables
uniform sampler2D uSampler;
uniform vec3 uKd;
uniform vec3 uKs;
uniform vec3 uLightPos;
uniform vec3 uCameraPos;
uniform vec3 uLightIntensity;

varying highp vec2 vTextureCoord;
varying highp vec3 vFragPos;
varying highp vec3 vNormal;

// Shadow map related variables
#define NUM_SAMPLES 40
#define BLOCKER_SEARCH_NUM_SAMPLES 40
#define PCF_NUM_SAMPLES 40
#define NUM_RINGS 10

#define EPS 1e-3
#define PI 3.141592653589793
#define PI2 6.283185307179586

uniform sampler2D uShadowMap;

varying vec4 vPositionFromLight;

highp float rand_1to1(highp float x ) { 
  // -1 -1
  return fract(sin(x)*10000.0);
}

highp float rand_2to1(vec2 uv ) { 
  // 0 - 1
	const highp float a = 12.9898, b = 78.233, c = 43758.5453;
	highp float dt = dot( uv.xy, vec2( a,b ) ), sn = mod( dt, PI );
	return fract(sin(sn) * c);
}

float unpack(vec4 rgbaDepth) {
    const vec4 bitShift = vec4(1.0, 1.0/255.0, 1.0/(255.0*255.0), 1.0/(255.0*255.0*255.0));
    return dot(rgbaDepth, bitShift);
}

vec2 poissonDisk[NUM_SAMPLES];

void poissonDiskSamples( const in vec2 randomSeed ) {
  float ANGLE_STEP = 3.883222077450933;// (sqrt(5)-1)/2 *2PI
  //float ANGLE_STEP = PI2 * float( NUM_RINGS ) / float( NUM_SAMPLES );
  float INV_NUM_SAMPLES = 1.0 / float( NUM_SAMPLES );

  float angle = rand_2to1( randomSeed ) * PI2;
  float radius = INV_NUM_SAMPLES;
  float radiusStep = radius;

  for( int i = 0; i < NUM_SAMPLES; i ++ ) {
    poissonDisk[i] = vec2( cos( angle ), sin( angle ) ) * pow( radius, 0.75 );
    radius += radiusStep;
    angle += ANGLE_STEP;
  }
}

void uniformDiskSamples( const in vec2 randomSeed ) {

  float randNum = rand_2to1(randomSeed);
  float sampleX = rand_1to1( randNum ) ;
  float sampleY = rand_1to1( sampleX ) ;

  float angle = sampleX * PI2;
  float radius = sqrt(sampleY);

  for( int i = 0; i < NUM_SAMPLES; i ++ ) {
    poissonDisk[i] = vec2( radius * cos(angle) , radius * sin(angle)  );

    sampleX = rand_1to1( sampleY ) ;
    sampleY = rand_1to1( sampleX ) ;

    angle = sampleX * PI2;
    radius = sqrt(sampleY);
  }
}

float PCF(sampler2D shadowMap, vec4 coords) {
  poissonDiskSamples(coords.xy);
  const float bias = 0.005;
  float sum = 0.0;
  for(int i = 0;i<PCF_NUM_SAMPLES;++i){
    float depthInShadowmap = unpack(texture2D(shadowMap,coords.xy+poissonDisk[i]*0.001).rgba);
    sum += ((depthInShadowmap + bias)< coords.z?0.0:1.0);
  }
  return sum/float(PCF_NUM_SAMPLES);
}


float findBlocker( sampler2D shadowMap, vec2 uv, float zReceiver ) {
  float dBlocker = zReceiver*0.01;
  const float wLight = 0.006;
  const float c = 10.0;
  float sampleSize = wLight*zReceiver*c;
  float sum = 0.01;
  for(int i = 0;i<BLOCKER_SEARCH_NUM_SAMPLES;++i){
    float depthInShadowmap = unpack(texture2D(shadowMap,uv+poissonDisk[i]*sampleSize).rgba);
    if(depthInShadowmap < zReceiver){
      dBlocker += depthInShadowmap;
      sum += 1.0;
    }
  }
  return dBlocker/float(sum);
}

float PCSS(sampler2D shadowMap, vec4 coords){
  poissonDiskSamples(coords.xy);
  // STEP 1: avgblocker depth
  float dBlocker = findBlocker(shadowMap,coords.xy,coords.z);
  // STEP 2: penumbra size
  const float wLight = 0.006;
  float penumbra = (coords.z-dBlocker)/dBlocker * wLight;
  // STEP 3: filtering
  const float bias = 0.005;
  float sum = 0.0;
  for(int i = 0;i<PCF_NUM_SAMPLES;++i){
    float depthInShadowmap = unpack(texture2D(shadowMap,coords.xy+poissonDisk[i]*penumbra).rgba);
    sum += ((depthInShadowmap + bias)< coords.z?0.0:1.0);
  }
  return sum/float(PCF_NUM_SAMPLES);
}


float useShadowMap(sampler2D shadowMap, vec4 shadowCoord){
  float depthInShadowmap = unpack(texture2D(shadowMap,shadowCoord.xy).rgba);
  const float bias = 0.005;
  if(depthInShadowmap + bias <shadowCoord.z)
    return 0.0;
  return 1.0;
}

vec3 blinnPhong() {
  vec3 color = texture2D(uSampler, vTextureCoord).rgb;
  color = pow(color, vec3(2.2));

  vec3 ambient = 0.05 * color;

  vec3 lightDir = normalize(uLightPos);
  vec3 normal = normalize(vNormal);
  float diff = max(dot(lightDir, normal), 0.0);
  vec3 light_atten_coff =
      uLightIntensity / pow(length(uLightPos - vFragPos), 2.0);
  vec3 diffuse = diff * light_atten_coff * color;

  vec3 viewDir = normalize(uCameraPos - vFragPos);
  vec3 halfDir = normalize((lightDir + viewDir));
  float spec = pow(max(dot(halfDir, normal), 0.0), 32.0);
  vec3 specular = uKs * light_atten_coff * spec;

  vec3 radiance = (ambient + diffuse + specular);
  vec3 phongColor = pow(radiance, vec3(1.0 / 2.2));
  return phongColor;
}

void main(void) {
  float visibility;
  vec3 projCoords = vPositionFromLight.xyz / vPositionFromLight.w;
  vec3 shadowCoord = projCoords * 0.5 + 0.5;
  //visibility = useShadowMap(uShadowMap, vec4(shadowCoord, 1.0));
  //visibility = PCF(uShadowMap, vec4(shadowCoord, 1.0));
  visibility = PCSS(uShadowMap, vec4(shadowCoord, 1.0));

  vec3 phongColor = blinnPhong();

  gl_FragColor = vec4(phongColor * visibility, 1.0);
}