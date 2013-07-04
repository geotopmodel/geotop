CC = gcc
CFLAGS = -I../fluidturtle -I ../../geotop

OBJS	= ${SRCS:.c=.o}

all: $(OBJLIB)

.c.o:
	$(CC) -c $(CFLAGS) $(CPPFLAGS) -o $@ $<

$(OBJLIB): $(OBJS)
	$(ECHO) $(AR) rcv $(OBJLIB) $?
	$(AR) rcv $(OBJLIB) $?

clean :
	-$(RM) -f $(OBJS) $(OBJLIB)
