# For Linux or any machines with gcc compiler
CC = gcc 
CFLAGS = -std=c99 -Wall -pedantic
BIN = ../bin

all: ImageReadWriteExample SurrogateFunctionExample SolveExample NonCausalPredictionError ICD_optimization ICD_optimization_blurred_and_noisy

clean:
	/bin/rm *.o $(BIN)/*

OBJ = tiff.o allocate.o randlib.o qGGMRF.o solve.o

ImageReadWriteExample: ImageReadWriteExample.o $(OBJ) 
	$(CC) $(CFLAGS) -o ImageReadWriteExample ImageReadWriteExample.o $(OBJ) -lm
	mv ImageReadWriteExample $(BIN)

SurrogateFunctionExample: SurrogateFunctionExample.o $(OBJ) 
	$(CC) $(CFLAGS) -o SurrogateFunctionExample SurrogateFunctionExample.o $(OBJ) -lm
	mv SurrogateFunctionExample $(BIN)

SolveExample: SolveExample.o $(OBJ) 
	$(CC) $(CFLAGS) -o SolveExample SolveExample.o $(OBJ) -lm
	mv SolveExample $(BIN)

NonCausalPredictionError: non_causal_prediction_error.o $(OBJ)
	$(CC) $(CFLAGS) -o ncpe non_causal_prediction_error.o $(OBJ) -lm
	mv ncpe $(BIN)

ICD_optimization: ICD_optimization.o $(OBJ)
	$(CC) $(CFLAGS) -o ICD_opt ICD_optimization.o $(OBJ) -lm
	mv ICD_opt $(BIN)

ICD_optimization_blurred_and_noisy: ICD_optimization_blurred_and_noisy.o $(OBJ)
	$(CC) $(CFLAGS) -o ICD_BN_opt ICD_optimization_blurred_and_noisy.o $(OBJ) -lm
	mv ICD_BN_opt $(BIN)

