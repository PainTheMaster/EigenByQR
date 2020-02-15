package main

import (
	"fmt"
	"math"
)

func main() {

	A := [][]float64{
		{16, -1, 1, 2},
		{2, 12, 1, -1},
		{1, 3, -24, 2},
		{4, -2, 1, 20}}

	y := []float64{20, 13, -16, 43}
	/*
		_ = eigenByQR(A, 1000)
	*/
	/*	x := solver(Qt, R, y)

		fmt.Println("x:\n", x)*/
	/*
		LU := lu(A)
		fmt.Println("LU:\n", LU)*/

	x := luSolver(A, y)
	fmt.Println("x:", x)
}

func eigenByQR(A [][]float64, rep int) (R [][]float64) {
	size := len(A)

	R = make([][]float64, size)
	for i := 0; i <= size-1; i++ {
		R[i] = make([]float64, size)
		for j := 0; j <= size-1; j++ {
			R[i][j] = A[i][j]
		}
	}

	fmt.Println("A:\n", A)
	fmt.Println("R:\n", R)

	var Qt [][]float64

	for i := 1; i <= rep; i++ {
		/*まずQt,Rを求める*/
		fmt.Println("pre-qr R:\n", R)

		Qt, R = qr(R) /*宣言同時代入(:=)してはいけない。そうしてしまう右辺と左辺のRが別々のシンボルとなり、毎回Rの初期値から初めてしまうようだ!*/

		/* Qを順に掛けていく */
		lenQt := len(Qt)
		for idxQt := 0; idxQt <= lenQt-1; idxQt++ {
			for row := 0; row <= size-1; row++ {
				helperMultiplyEigenHorizon(R, row, Qt[idxQt])
			}
		}

		fmt.Println(i, ", R:\n", R)
	}

	return
}

func helperMultiplyEigenHorizon(R [][]float64, row int, vt []float64) {

	sizeMat := len(R)
	sizeVec := len(vt)

	var innerprod float64

	innerprod = 0.0
	for j := 0; j <= sizeVec-1; j++ {
		innerprod += R[row][sizeMat-sizeVec+j] * vt[j]
	}

	for j := 0; j <= sizeVec-1; j++ {
		R[row][sizeMat-sizeVec+j] = R[row][sizeMat-sizeVec+j] - 2.0*vt[j]*innerprod
	}
}

func qr(A [][]float64) (Qt [][]float64, R [][]float64) {

	size := len(A)

	Qt = make([][]float64, size-1)

	R = make([][]float64, size)
	for i := 0; i <= size-1; i++ {
		R[i] = make([]float64, size)
	}
	for i := 0; i <= size-1; i++ {
		for j := 0; j <= size-1; j++ {
			R[i][j] = A[i][j]
		}
	}

	for colPiv := 0; colPiv <= size-2; colPiv++ {
		tempH, topVal := householder(R, colPiv)

		Qt[colPiv] = make([]float64, size-colPiv)
		for row := 0; row <= size-colPiv-1; row++ {
			Qt[colPiv][row] = tempH[row]
		}

		R[colPiv][colPiv] = topVal
		for i := colPiv + 1; i <= size-1; i++ {
			R[i][colPiv] = 0.0
		}

		for colProduct := colPiv + 1; colProduct <= size-1; colProduct++ {
			helperMultiplyQrVertical(R, colProduct, tempH)
		}

	}

	return
}

func helperMultiplyQrVertical(R [][]float64, col int, v []float64) {
	sizeMat := len(R)
	sizeVec := len(v)

	var innerprod float64
	rowTop := sizeMat - sizeVec

	innerprod = 0.0
	for i := 0; i <= sizeVec-1; i++ {
		innerprod += R[rowTop+i][col] * v[i]
	}

	for i := 0; i <= sizeVec-1; i++ {
		R[rowTop+i][col] = R[rowTop+i][col] - 2.0*v[i]*innerprod
	}
}

//Gives Householder "vector" of A focusing on the column c "col" and topVal
func householder(A [][]float64, col int) (h []float64, topVal float64) {
	size := len(A)

	topVal = 0.0
	for i := col; i <= size-1; i++ {
		topVal += A[i][col] * A[i][col]
	}

	tempSqNormH := topVal

	topVal = math.Sqrt(topVal)
	if A[col][col] > 0 {
		topVal *= -1.0
	}

	h = make([]float64, size-col)

	h[0] = A[col][col] - topVal

	for i := 1; i <= size-col-1; i++ {
		h[i] = A[col+i][col]
	}

	tempSqNormH -= A[col][col] * A[col][col]
	tempSqNormH += h[0] * h[0]

	normalizFactor := 1.0 / math.Sqrt(tempSqNormH)

	for i := 0; i <= size-col-1; i++ {
		h[i] *= normalizFactor
	}

	return
}

func qrSolver(Qt [][]float64, R [][]float64, y []float64) (x []float64) {

	size := len(R)

	for idxQ := 0; idxQ <= size-2; idxQ++ {
		innerProd := 0.0
		for i := 0; i+idxQ <= size-1; i++ {
			innerProd += Qt[idxQ][i] * y[idxQ+i]
		}

		for i := 0; i+idxQ <= size-1; i++ {
			y[idxQ+i] = y[idxQ+i] - 2.0*Qt[idxQ][i]*innerProd
		}
	}

	x = make([]float64, size)

	for rowSolv := size - 1; rowSolv >= 0; rowSolv-- {
		var sum float64
		for i := size - 1; i >= rowSolv+1; i-- {
			sum += R[rowSolv][i] * x[i]
		}
		x[rowSolv] = (y[rowSolv] - sum) / R[rowSolv][rowSolv]
	}

	return

}

//Lの対角成分が1であるとする。
func lu(A [][]float64) (LU [][]float64) {
	size := len(A)

	LU = make([][]float64, size)
	for i := 0; i <= size-1; i++ {
		LU[i] = make([]float64, size)
	}

	for col := 0; col <= size-1; col++ {
		for row := 0; row <= col; row++ {
			var innerprod float64
			for i := 0; i <= row-1; i++ {
				innerprod += LU[row][i] * LU[i][col]
			}
			LU[row][col] = A[row][col] - innerprod
		}

		for row := col + 1; row <= size-1; row++ {
			var innerprod float64
			for i := 0; i <= row-1; i++ {
				innerprod += LU[row][i] * LU[i][col]
			}
			LU[row][col] = (A[row][col] - innerprod) / LU[col][col]
		}
	}
	return
}

//改定必必要。Lの対角成分を1にしなければならない
func luInverse(LU [][]float64) (LUinv [][]float64) {

	size := len(LU)

	LUinv = make([][]float64, size)
	for i := 0; i <= size-1; i++ {
		LUinv[i] = make([]float64, size)
	}

	//U
	for col := size - 1; col >= 0; col-- {
		row := col
		LUinv[col][col] = 1.0 / LU[col][col]
		row--
		for ; row >= 0; row-- {
			var innerprod float64
			for i := row + 1; i <= col; i++ {
				innerprod += LU[row][i] * LUinv[i][col]
			}
			//			innerprod += LU[row][col] * 1.0
			LUinv[row][col] = -1.0 * innerprod / LU[row][row]
		}
	}

	//L
	for row := size - 1; row >= 0; row-- {

		for col := row - 1; col >= 0; col-- {
			var innerprod float64
			for i := col + 1; i <= row-1; i++ {
				innerprod += LUinv[row][i] * LU[i][col]
			}
			innerprod += 1.0 * LU[row][col]
			LUinv[row][col] = -1.0 * innerprod
		}
	}
	return
}

func inverseMatrix(A [][]float64) (AInv [][]float64) {
	size := len(A)
	AInv = make([][]float64, size)
	for i := 0; i <= size-1; i++ {
		AInv[i] = make([]float64, size)
	}
	LU := lu(A)
	LUinv := luInverse(LU)

	//UL
	for row := 0; row <= size-1; row++ {
		//colが小さい時（長い列ベクトルでUが制限因子になっている時）
		for col := 0; col <= row-1; col++ {
			var innerprod float64
			for i := row; i <= size-1; i++ {
				innerprod += LUinv[row][i] * LUinv[i][col]
			}
			AInv[row][col] = innerprod
		}

		for col := row; col <= size-1; col++ {
			var innerprod float64
			innerprod = LUinv[row][col] * 1.0
			for i := col + 1; i <= size-1; i++ {
				innerprod += LUinv[row][i] * LUinv[i][col]
			}
			AInv[row][col] = innerprod
		}
	}

	return
}

func luSolver(A [][]float64, y []float64) (x []float64) {
	size := len(A)
	x = make([]float64, size)
	interm := make([]float64, size)

	LU := lu(A)
	//L*interm = y, interm = Ux
	for i := 0; i <= size-1; i++ {
		var innerprod float64
		for j := 0; j <= i-1; j++ {
			innerprod += LU[i][j] * interm[j]
		}
		interm[i] = y[i] - innerprod
	}
	//Ux = interm
	for i := size - 1; i >= 0; i-- {
		var innerprod float64
		for j := size - 1; j >= i+1; j-- {
			innerprod += LU[i][j] * x[j]
		}
		x[i] = (interm[i] - innerprod) / LU[i][i]
	}
	return
}
