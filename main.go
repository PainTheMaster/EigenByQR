package main

import "math"

func main() {

}

func qr(A [][]float64, rep int) (Qt [][]float64, R [][]float64) {

	size := len(A)

	Qt = make([][]float64, size-1)
	for i := 0; i <= size-1; i++ {
		Qt[i] = make([]float64, size-i)
	}

	R = make([][]float64, size)
	for i := 0; i <= size-1; i++ {
		R[i] = make([]float64, size)
	}

	{
		colPiv := 0
		tempH, topVal := householder(A, colPiv)
		for row := 0; row <= size-colPiv; row++ {
			Qt[colPiv][row] = tempH[row]
		}

		R[colPiv][colPiv] = topVal

		vecTemp := make([]float64, size-colPiv)
		for colProd := colPiv + 1; colProd <= size-1; colProd++ {
			innerprod := 0.0
			for i := colPiv; i <= size-1; i++ {
				innerprod += tempH[i-colPiv] * A[i][colProd]
			}

			for rowProd := colPiv; rowProd <= size-1; rowProd++ {
				vecTemp[rowProd-colPiv] = A[rowProd][colProd] - 2.0*tempH[rowProd-colPiv]*innerprod
			}
			for rowCopy := colPiv; rowCopy <= size-1; rowCopy++ {
				R[rowCopy][colProd] = vecTemp[rowCopy-colPiv]
			}
		}
	}

	for colPiv := 1; colPiv <= size-2; colPiv++ {
		tempH, topVal := householder(A, colPiv)
		for row := 0; row <= size-colPiv; row++ {
			Qt[colPiv][row] = tempH[row]
		}

		R[colPiv][colPiv] = topVal

		vecTemp := make([]float64, size-colPiv)
		for colProd := colPiv + 1; colProd <= size-1; colProd++ {
			innerprod := 0.0
			for i := colPiv; i <= size-1; i++ {
				innerprod += tempH[i-colPiv] * R[i][colProd]
			}

			for rowProd := colPiv; rowProd <= size-1; rowProd++ {
				vecTemp[rowProd-colPiv] = R[rowProd][colProd] - 2.0*tempH[rowProd-colPiv]*innerprod
			}
			for rowCopy := colPiv; rowCopy <= size-1; rowCopy++ {
				R[rowCopy][colProd] = vecTemp[rowCopy-colPiv]
			}
		}
	}

	return

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

	h[0] = A[0][col] - topVal

	for i := 1; i <= size-col-1; i++ {
		h[i] = A[i][col]
	}

	tempSqNormH -= A[col][col] * A[col][col]
	tempSqNormH += h[0] * h[0]

	normalizFactor := 1.0 / math.Sqrt(tempSqNormH)

	for i := 0; i <= size-col-1; i++ {
		h[i] /= normalizFactor
	}

	return

}
