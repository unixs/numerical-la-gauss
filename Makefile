PROG = gauss

CFLAGS = -Wall -Wextra -Wpedantic -ggdb -std=gnu11 -O0 `pkg-config --cflags glib-2.0`
LDLIBS = `pkg-config --libs glib-2.0 gsl`

default: bin/$(PROG)

bin/$(PROG): bin/$(PROG).o
	$(CC) $(CFLAGS) $^ $(LDLIBS) -o $@

bin/$(PROG).o: src/$(PROG).c
	$(CC) $(CFLAGS) -o $@ -c $<

clean:
	rm -fv bin/*.o bin/$(PROG)
