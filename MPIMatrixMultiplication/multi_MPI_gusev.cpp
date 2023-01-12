#include <clocale>
#include <stdlib.h>
#define _CRTDBG_MAP_ALLOC
#include <windows.h>
#include <crtdbg.h>
#include "mpi.h"
#include <iostream>

bool Test = false;

class manipulationMatrix
{
public:
	int* _matrix = NULL;

public:
	int row = 0;
	int column = 0; 


public:
	manipulationMatrix(int rows, int cols, bool empty = true)
	{
		row = rows;
		column = cols;
		generateMatrix(rows, cols, empty);
	}

	~manipulationMatrix()
	{
		// Очистка
		if (_matrix != NULL) delete[] _matrix;
		_matrix = NULL;
	}
private:
	// Создание матрицы со случайными числами
	void generateMatrix(int rows, int cols, bool empty)
	{
		_matrix = new int[rows * cols];

		for (int row = 0; row < rows; row++)
		{
			for (int col = 0; col < cols; col++)
			{
				int initNumber = 0;
				if (empty == false) initNumber = 1 + std::rand() % 10;

				_matrix[(row * cols + col)] = initNumber;
			}
		}
	}

public:
	void setElement(int row, int col, int value)
	{
		_matrix[(row * column + col)] = value;
	}

	int getElement(int row, int col)
	{
		return _matrix[(row * column + col)];
	}

	void printStr(int* buffer, int countrow, int countcols)
	{
		for (int r = 0; r < countrow; r++)
		{
			for (int j = 0; j < countcols; j++)
			{
				int index = r * countcols + j;
				std::cout << buffer[index] << " ";
			}
			std::cout << "\n";
		}
	}

	void printMatrix()
	{
		for (int r = 0; r < row; r++)
		{
			for (int c = 0; c < column; c++)
			{
				std::cout.width(4);
				std::cout << getElement(r, c);
			}
			std::cout << "\n";
		}
	}
};
int M = 2000;
int N = 1500;

void multiplactionSegmentA(int* rowsMatrixA, int* outbuffer, int countrowsMatrixA, manipulationMatrix* B, manipulationMatrix* matrixC)
{
	for (int mArow = 0; mArow < countrowsMatrixA; mArow++)
	{
		for (int Bcol = 0; Bcol < M; Bcol++)
		{
			int C = 0;
			for (int BrowAcol = 0; BrowAcol < N; BrowAcol++)
			{
				int t = mArow * N + BrowAcol;
				C += rowsMatrixA[t] * (*B).getElement(BrowAcol, Bcol);
			}
			if (outbuffer != NULL) outbuffer[mArow * M + Bcol] = C;
			if (matrixC != NULL) (*matrixC).setElement(mArow, Bcol, C);
		}
	}
}


void GetRows(int numrowsMatrixA, int rank, int worldsize, int* indexbegin, int* countrows)
{
	int d = numrowsMatrixA % worldsize;
	int wr = (numrowsMatrixA - d) / worldsize;
	int bgn = 0;
	if (rank > 0) bgn = wr * rank + d;
	int rows = wr + d;
	if (rank > 0) rows = wr;
	*indexbegin = bgn;
	*countrows = rows;
}


// Создание матриц, вычисление матрицы С, просмотр рез.
void WorkRoot(int world_rank, int world_size)
{
	std::cout << "Перемножение матриц с помощью mpi\n";
	std::cout << "Начнем подсчет!\n";
	double startInit = MPI_Wtime();
	std::srand(std::time(nullptr));
	manipulationMatrix matrixA = manipulationMatrix(M, N, false);
	manipulationMatrix matrixB = manipulationMatrix(N, M, false);
	manipulationMatrix matrixC(M, M);

	std::cout << "\nМатрица A: строки=" << matrixA.row << ", колонки=" << matrixA.column << "\n";
	std::cout << "Матрица B: строки=" << matrixB.row << ", колонки=" << matrixB.column << "\n";
	std::cout << "Матрица C: строки=" << matrixC.row << ", колонки=" << matrixC.column << "\n";
	std::cout << "\n";

	double endInit = MPI_Wtime();
	double t1 = endInit - startInit;
	double startSend = MPI_Wtime();
	for (int i = 1; i < world_size; i++)
	{
		int id = 0, rowCount = 0;
		GetRows(M, i, world_size, &id, &rowCount);
		int countElements = rowCount * N;
		int offsetBuffer = id * N;
		int error1 = MPI_Send(matrixA._matrix + offsetBuffer, countElements, MPI_INT, i, 0, MPI_COMM_WORLD);
		int error = MPI_Send(matrixB._matrix, N * M, MPI_INT, i, 0, MPI_COMM_WORLD);
	}

	double endSend = MPI_Wtime();
	double timeSent = endSend - startSend; 
	// Умножение соответствующего количества строк матрицы А и столбцов матрицы В
	std::cout << "\nПеремножение (rank 0)\n";
	double startMultiply = MPI_Wtime();
	int index = 0, countrows = 0;
	GetRows(M, world_rank, world_size, &index, &countrows);

	std::cout << "rank " << world_rank << " кол. строк = " << countrows << "\n";

	// Умножение матриц строками-сегментами матрицы А
	multiplactionSegmentA(matrixA._matrix + index * N, NULL, countrows, &matrixB, &matrixC);
	double endMultiply = MPI_Wtime();
	double timeMultiplication = endMultiply - startMultiply;
	double startRecv = MPI_Wtime();

	for (int i = 1; i < world_size; i++)
	{
		int index = 0, countrows = 0;
		int rank = i;
		int size_comm = world_size;

		GetRows(M, rank, size_comm, &index, &countrows);

		int* recvRowsMatrixC = new int[countrows * M];
		memset(recvRowsMatrixC, 0, countrows * M);

		int errorRecv = MPI_Recv(recvRowsMatrixC, countrows * M, MPI_INT, i, 0, MPI_COMM_WORLD,
			MPI_STATUS_IGNORE);
		memcpy(matrixC._matrix + (index * M), recvRowsMatrixC, countrows * M * sizeof(int));

		delete[] recvRowsMatrixC;
	}

	double endRecv = MPI_Wtime();
	double timeReceiveResult = endRecv - startRecv; 
	std::cout << "\n\n";
	std::cout << "Итоговое время:" << t1 + timeSent + timeMultiplication + timeReceiveResult;
}


void processorJob(int world_rank, int world_size)
{
	double startRecvMatrixSource = MPI_Wtime();
	manipulationMatrix recvB = manipulationMatrix(N, M);

	int indexBeginRows = 0;
	int countrowsMatrixA = 0;
	GetRows(M, world_rank, world_size, &indexBeginRows, &countrowsMatrixA);
	int* recvRowsMatrixA = new int[countrowsMatrixA * N];

	int error1 = MPI_Recv(recvRowsMatrixA, countrowsMatrixA * N, MPI_INT, 0, 0, MPI_COMM_WORLD,
		MPI_STATUS_IGNORE);


	int error = MPI_Recv(recvB._matrix, N * M, MPI_INT, 0, 0, MPI_COMM_WORLD,
		MPI_STATUS_IGNORE);

	double endRecvMatrixSource = MPI_Wtime();
	double timeResultRecvMatrixSource = endRecvMatrixSource - startRecvMatrixSource; // 1
	std::cout << "Time recv: " << timeResultRecvMatrixSource << " sec.\n";
	double startMultiply = MPI_Wtime();
	int* outRowsMatrixC = new int[countrowsMatrixA * M];
	multiplactionSegmentA(recvRowsMatrixA, outRowsMatrixC, countrowsMatrixA, &recvB, NULL);
	delete[] recvRowsMatrixA;
	double endMultiply = MPI_Wtime();
	double timeResultMultiply = endMultiply - startMultiply; // 2
	std::cout << "Time multiply: " << timeResultMultiply << " sec.\n";
	double startSend = MPI_Wtime();

	// Отправка результата для матрицы С
	int errorSend = MPI_Send(outRowsMatrixC, countrowsMatrixA * M, MPI_INT, 0, 0, MPI_COMM_WORLD);
	delete[] outRowsMatrixC;
	double endSend = MPI_Wtime();
	double timeResultSendRows = endSend - startSend; 
	std::cout << "Time send: " << timeResultSendRows << " sec.\n";
	std::cout << "Итог: " << timeResultRecvMatrixSource + timeResultMultiply + timeResultSendRows;
}



int main(int argc, char* argv[])
{
	setlocale(LC_ALL, "Russian");

	if (argc > 1 &&
		(strcmp(argv[1], "Test") == 0 || strcmp(argv[1], "test") == 0))
	{
		Test = true;
	}
	MPI_Init(NULL, NULL);
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
	if (world_rank == 0)
	{
		WorkRoot(world_rank, world_size);
	}
	else
	{
		processorJob(world_rank, world_size);
	}
	MPI_Finalize();
	_CrtDumpMemoryLeaks();
}
