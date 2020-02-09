package main

import (
	"fmt"
	"math"
)

func main() {

	A := [][]float64{
		{1, 1, 1},
		{1, 2, 4},
		{1, 3, 6}}

	y := []float64{3, 7, 10}

	Qt, R := qr(A)

	x := solver(Qt, R, y)

	fmt.Println("R:\n", R)

	fmt.Println()
	fmt.Println("x:\n", x)
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

func solver(Qt [][]float64, R [][]float64, y []float64) (x []float64) {

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
