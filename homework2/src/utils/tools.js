function getRotationPrecomputeL(precompute_L, rotationMatrix){
	let rotMat = mat4Matrix2mathMatrix(rotationMatrix);
	let rotateSHMat1 = computeSquareMatrix_3by3(rotMat);
	let rotateSHMat2 = computeSquareMatrix_5by5(rotMat);
	let ret = [];
	for(i = 0; i < 3; i++)
	{
		let colors = math.clone(precompute_L[i]);
		let sh1 = [colors[1],colors[2],colors[3]];
		let sh2 = [colors[4],colors[5],colors[6],colors[7],colors[8]];
		let rotatedSH1 = math.multiply(rotateSHMat1, [sh1[0], sh1[1], sh1[2]]);
		let rotatedSH2 = math.multiply(rotateSHMat2, [sh2[0], sh2[1], sh2[2], sh2[3], sh2[4]]);
		colors[1] = rotatedSH1[0];
		colors[2] = rotatedSH1[1];
		colors[3] = rotatedSH1[2];
		colors[4] = rotatedSH2[0];
		colors[5] = rotatedSH2[1];
		colors[6] = rotatedSH2[2];
		colors[7] = rotatedSH2[3];
		colors[8] = rotatedSH2[4];
		ret.push([colors[0], colors[1], colors[2], 
				  colors[3], colors[4], colors[5],
				  colors[6], colors[7], colors[8]]);
	}
	return ret;
}

function computeSquareMatrix_3by3(rotationMatrix){ 
	// 计算方阵SA(-1) 3*3 
	// 1、pick ni - {ni}
	let n1 = [1, 0, 0, 0]; let n2 = [0, 0, 1, 0]; let n3 = [0, 1, 0, 0];
	// 2、{P(ni)} - A
	let Pn1 = SHEval(n1[0],n1[1],n1[2],3);
	let Pn2 = SHEval(n2[0],n2[1],n2[2],3);
	let Pn3 = SHEval(n3[0],n3[1],n3[2],3);
	let A = math.matrix([
		[Pn1[1],Pn2[1],Pn3[1]],
		[Pn1[2],Pn2[2],Pn3[2]],
		[Pn1[3],Pn2[3],Pn3[3]]
	]);
	// 3、用 R 旋转 ni - {R(ni)}
	let Rn1 = math.multiply(rotationMatrix,n1);
	let Rn2 = math.multiply(rotationMatrix,n2);
	let Rn3 = math.multiply(rotationMatrix,n3);
	// 4、R(ni) SH投影 - S
	let PRn1 = SHEval(Rn1[0], Rn1[1], Rn1[2], 3);
	let PRn2 = SHEval(Rn2[0], Rn2[1], Rn2[2], 3);
	let PRn3 = SHEval(Rn3[0], Rn3[1], Rn3[2], 3);
	let S = math.matrix([
		[PRn1[1],PRn2[1],PRn3[1]],
		[PRn1[2],PRn2[2],PRn3[2]],
		[PRn1[3],PRn2[3],PRn3[3]]
	]);
	// 5、S*A_inverse
	return math.transpose(math.multiply(S._data,math.inv(A)._data));
}

function computeSquareMatrix_5by5(rotationMatrix){ 
	// 计算方阵SA(-1) 5*5
	// 1、pick ni - {ni}	
	let k = 1 / math.sqrt(2);
	let n1 = [1, 0, 0, 0]; let n2 = [0, 0, 1, 0]; let n3 = [k, k, 0, 0]; 
	let n4 = [k, 0, k, 0]; let n5 = [0, k, k, 0];
	// 2、{P(ni)} - A
	let Pn1 = SHEval(n1[0],n1[1],n1[2],5);
	let Pn2 = SHEval(n2[0],n2[1],n2[2],5);
	let Pn3 = SHEval(n3[0],n3[1],n3[2],5);
	let Pn4 = SHEval(n4[0],n4[1],n4[2],5);
	let Pn5 = SHEval(n5[0],n5[1],n5[2],5);	
	let A = math.matrix([
		[Pn1[4],Pn2[4],Pn3[4],Pn4[4],Pn5[4]],
		[Pn1[5],Pn2[5],Pn3[5],Pn4[5],Pn5[5]],
		[Pn1[6],Pn2[6],Pn3[6],Pn4[6],Pn5[6]],
		[Pn1[7],Pn2[7],Pn3[7],Pn4[7],Pn5[7]],
		[Pn1[8],Pn2[8],Pn3[8],Pn4[8],Pn5[8]]
	]);
	// 3、用 R 旋转 ni - {R(ni)}
	let Rn1 = math.multiply(rotationMatrix,n1);
	let Rn2 = math.multiply(rotationMatrix,n2);
	let Rn3 = math.multiply(rotationMatrix,n3);
	let Rn4 = math.multiply(rotationMatrix,n4);
	let Rn5 = math.multiply(rotationMatrix,n5);
	// 4、R(ni) SH投影 - S
	let PRn1 = SHEval(Rn1[0], Rn1[1], Rn1[2], 5);
	let PRn2 = SHEval(Rn2[0], Rn2[1], Rn2[2], 5);
	let PRn3 = SHEval(Rn3[0], Rn3[1], Rn3[2], 5);
	let PRn4 = SHEval(Rn4[0], Rn4[1], Rn4[2], 5);
	let PRn5 = SHEval(Rn5[0], Rn5[1], Rn5[2], 5);
	let S = math.matrix([
		[PRn1[4],PRn2[4],PRn3[4],PRn4[4],PRn5[4]],
		[PRn1[5],PRn2[5],PRn3[5],PRn4[5],PRn5[5]],
		[PRn1[6],PRn2[6],PRn3[6],PRn4[6],PRn5[6]],
		[PRn1[7],PRn2[7],PRn3[7],PRn4[7],PRn5[7]],
		[PRn1[8],PRn2[8],PRn3[8],PRn4[8],PRn5[8]]
	]);
	// 5、S*A_inverse
	return math.transpose(math.multiply(S._data,math.inv(A)._data));
}

function mat4Matrix2mathMatrix(rotationMatrix){

	let mathMatrix = [];
	for(let i = 0; i < 4; i++){
		let r = [];
		for(let j = 0; j < 4; j++){
			r.push(rotationMatrix[i*4+j]);
		}
		mathMatrix.push(r);
	}
	//return math.matrix(mathMatrix)
	return mathMatrix
}

function getMat3ValueFromRGB(precomputeL){

    let colorMat3 = [];
    for(var i = 0; i<3; i++){
        colorMat3[i] = mat3.fromValues( precomputeL[0][i], precomputeL[1][i], precomputeL[2][i],
										precomputeL[3][i], precomputeL[4][i], precomputeL[5][i],
										precomputeL[6][i], precomputeL[7][i], precomputeL[8][i] ); 
	}
    return colorMat3;
}