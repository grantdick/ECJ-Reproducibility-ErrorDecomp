CC:=gcc
SRCDIR:=src
OBJDIR:=build
INCDIR:=include
BINDIR:=dist

INCS:=$(wildcard $(SRCDIR)/*.h)
SRCS:=$(wildcard $(SRCDIR)/*.c)
OBJS:=$(subst $(SRCDIR)/,$(OBJDIR)/,$(patsubst %.c,%.o,$(wildcard $(SRCDIR)/*.c)))

COMMON_OBJS:=$(OBJDIR)/baselearner.o $(OBJDIR)/data.o \
	$(OBJDIR)/parsetree.o $(OBJDIR)/problemgenerator.o \
	$(OBJDIR)/readline.o $(OBJDIR)/rng.o $(OBJDIR)/util.o

CFLAGS:=-std=gnu89 -Wall -Wextra -pedantic -march=native -O3 -g
IFLAGS:=-I$(INCDIR)
LFLAGS:=-lm

BAGGP_BIN:=$(BINDIR)/baggp
GSGP_BIN:=$(BINDIR)/gsgp
SSPBE_BIN:=$(BINDIR)/sspbe
SSEGP_BIN:=$(BINDIR)/2segp

all: $(BAGGP_BIN) $(GSGP_BIN) $(SSPBE_BIN) $(SSEGP_BIN)

$(BAGGP_BIN): $(COMMON_OBJS) $(OBJDIR)/baggp.o $(OBJDIR)/baggp_main.o
	@echo linking $@ from $^
	@mkdir -p $(BINDIR)
	@$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

$(GSGP_BIN): $(COMMON_OBJS) $(OBJDIR)/gsgp.o $(OBJDIR)/gsgp_main.o
	@echo linking $@ from $^
	@mkdir -p $(BINDIR)
	@$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

$(SSPBE_BIN): $(COMMON_OBJS) $(OBJDIR)/sspbe.o $(OBJDIR)/sspbe_ensemble.o $(OBJDIR)/sspbe_main.o
	@echo linking $@ from $^
	@mkdir -p $(BINDIR)
	@$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

$(SSEGP_BIN): $(COMMON_OBJS) $(OBJDIR)/2segp.o $(OBJDIR)/2segp_ensemble.o $(OBJDIR)/2segp_main.o
	@echo linking $@ from $^
	@mkdir -p $(BINDIR)
	@$(CC) $(CFLAGS) $^ -o $@ $(LFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.c $(INCS)
	@echo compiling $< into $@
	@mkdir -p $(OBJDIR)
	@$(CC) $(CFLAGS) $(IFLAGS) -c $< -o $@

clean:
	@rm -rf $(OBJDIR)

nuke: clean
	@rm -rf $(INCDIR) $(BINDIR)

strip: all
	@echo running strip on $(BIN)
	@strip $(BIN)
