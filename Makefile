# Compiler
CC      = gcc

# Compiler flags
CFLAGS  = -Wall -Wextra -O3 -march=native -Iinclude

# Linker flags
LDFLAGS =
LIBS    = -lfftw3 -lm

# Directories
SRC_DIR = codes/src
OBJ_DIR = obj
BIN_DIR = exec

# Source files
SRC = $(wildcard $(SRC_DIR)/*.c)

# Object files
OBJ = $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRC))

# Executable
TARGET = $(BIN_DIR)/code1Dw.out

# ----------------------------

all: $(TARGET)

$(TARGET): $(OBJ)
	mkdir -p $(BIN_DIR)
	$(CC) $(OBJ) $(LDFLAGS) $(LIBS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OBJ_DIR) $(BIN_DIR)

run: $(TARGET)
	./$(TARGET)

.PHONY: all clean run
