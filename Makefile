TARGET = NL

CC = gcc
SOURCE_DIR = src
INCLUDE_DIR = $(SOURCE_DIR)/*/
CFLAGS = -I$(INCLUDE_DIR) -I$(SOURCE_DIR)/
C_SOURCES = $(wildcard src/*/*.c src/*.c *.c )
DEPS = $(wildcard src/*/.h src/*.h *.h)
OBJ = ${C_SOURCES:.c=.o}

# First rule is the one executed when no parameters are fed to the Makefile


%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(TARGET): $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)

run: $(TARGET)
	./$(TARGET)

clean:
	$(RM) $(INCLUDE_DIR)*.o $(INCLUDE_DIR)*.dis $(INCLUDE_DIR)*.elf $(INCLUDE_DIR)*.bin
	$(RM) src/*.bin src/*.o src/*.dis src/*.elf
	$(RM) lib/*.o
