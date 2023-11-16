TARGET = my_prog

BINDIR = bin
OBJDIR = obj

SOURCES = src/tilt_solver.c \
main.c

INC1 = /usr/local/include
INC2 = src

LIB = /usr/local/lib

ifeq ($(DEBUG), 1)
	override CFLAGS:=-g $(CFLAGS)
	override LDFLAGS:=-g $(LDFLAGS)
endif


WARNING = -Wall -Wextra -Wpedantic
CFLAGS = -c
LDFLAGS = -lgsl -lgslcblas -lm

BUILDDIR = $(BINDIR)/$(OBJDIR)


OBJECTS = $(addprefix $(BUILDDIR)/,$(patsubst %.c,%.o,$(notdir $(SOURCES))))
vpath %.c $(sort $(dir $(SOURCES)))

all: $(BUILDDIR) $(BINDIR)/$(TARGET)
	@echo Compiling done

$(BINDIR)/$(TARGET): $(OBJECTS)
	@echo LD       $@
	@gcc -L$(LIB) $(OBJECTS) $(LDFLAGS) -o $@

$(BUILDDIR)/%.o: %.c
	@echo CC       $^
	@gcc $(CFLAGS) $(WARNING) -I$(INC1) -I$(INC2) $< -o $@ 

$(BUILDDIR):
	@mkdir -p $@


clean:
	@rm -f $(BINDIR)/my_prog
	@rm -f -r $(BINDIR)
	@rm -f *.txt
	@echo "Cleaning done"

print:
	@echo $(OBJECTS)
	@echo $(SOURCES)
	@echo $(BINDIR)/$(TARGET)