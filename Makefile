PRINTDEBUG := 0

TARGET := sw_sam_sort
BIN_TARGET := $(TARGET)

DIR_SRC := ./src
SLAVE_DIR_SRC := ./slave
DIR_OBJ := ./obj

PREFIX ?= /usr/local
BINDIR ?= $(PREFIX)/bin

INCLUDE_DIRS ?=
LIBRARY_DIRS ?=
LIBS := -lathread

SRC := $(wildcard ${DIR_SRC}/*.c)
OBJ := $(patsubst %.c,${DIR_OBJ}/%.o,$(notdir ${SRC}))

SLAVE_SRC := $(wildcard ${SLAVE_DIR_SRC}/*.c)
SLAVE_OBJ := $(patsubst %.c,${DIR_OBJ}/slave_%.o,$(notdir ${SLAVE_SRC}))

OBJ_ALL := $(OBJ) $(SLAVE_OBJ)

CC_HOST = swgcc
CC_SLAVE = swgcc

CFLAGS_COMMON := -g -O3 -w
# 这里可以按需增加头文件搜索路径
CFLAGS_HOST := $(CFLAGS_COMMON) $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) -I./slave
CFLAGS_SLAVE := $(CFLAGS_COMMON) $(foreach includedir,$(INCLUDE_DIRS),-I$(includedir)) -I./slave

LDFLAGS := $(foreach librarydir,$(LIBRARY_DIRS),-L$(librarydir)) $(LIBS)

all: $(BIN_TARGET)

$(BIN_TARGET): $(OBJ_ALL)
	$(CC_HOST) -mhybrid $^ -o $@ $(LDFLAGS)

${DIR_OBJ}/%.o:${DIR_SRC}/%.c
	@mkdir -p $(DIR_OBJ)
	$(CC_HOST) -mhost -c $< -o $@ $(CFLAGS_HOST)

${DIR_OBJ}/slave_%.o:${SLAVE_DIR_SRC}/%.c
	@mkdir -p $(DIR_OBJ)
	$(CC_SLAVE) -mslave -msimd -c $< -o $@ $(CFLAGS_SLAVE)

.PHONY: all clean install

clean:
	rm -rf $(DIR_OBJ)/*.o
	rm -f $(TARGET)

install: $(BIN_TARGET)
	install $(BIN_TARGET) $(BINDIR)/$(TARGET)
	@echo "Installed."


