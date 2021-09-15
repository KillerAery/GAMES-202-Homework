#ifdef GL_ES
precision mediump float;
#endif

// Phong related variables
uniform sampler2D uSampler;
uniform vec3 uKd;
uniform vec3 uKs;
uniform vec3 uLightPos;
uniform vec3 uCameraPos;
uniform vec3 uLightRadiance;
uniform mat3 uPrecomputeLR;
uniform mat3 uPrecomputeLG;
uniform mat3 uPrecomputeLB;

varying highp vec2 vTextureCoord;
varying highp vec3 vFragPos;
varying highp vec3 vNormal;
varying highp mat3 vPrecomputeLT;


void main(void) {
  //vec3 color = texture2D(uSampler, vTextureCoord).rgb;
  vec3 color = vec3(
      dot(uPrecomputeLR[0],vPrecomputeLT[0])+dot(uPrecomputeLR[1],vPrecomputeLT[1])+dot(uPrecomputeLR[2],vPrecomputeLT[2]),
      dot(uPrecomputeLG[0],vPrecomputeLT[0])+dot(uPrecomputeLG[1],vPrecomputeLT[1])+dot(uPrecomputeLG[2],vPrecomputeLT[2]),
      dot(uPrecomputeLB[0],vPrecomputeLT[0])+dot(uPrecomputeLB[1],vPrecomputeLT[1])+dot(uPrecomputeLB[2],vPrecomputeLT[2])
  );
  gl_FragColor = vec4(0.75*color, 1.0);
}
