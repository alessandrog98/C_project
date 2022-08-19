main_iplib: main_iplib.o ip_lib.o bmp.o ip_lib.h bmp.h
	gcc main_iplib.o ip_lib.o bmp.o -omain_iplib -Wall -lm --ansi --pedantic -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

main_iplib.o: main_iplib.c ip_lib.h bmp.h
	gcc -c main_iplib.c -omain_iplib.o -Wall -lm --ansi --pedantic -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

ip_lib.o: ip_lib.c ip_lib.h bmp.h
	gcc -c ip_lib.c -oip_lib.o -Wall -lm --ansi --pedantic -g3 -O3 -fsanitize=address -fsanitize=undefined -std=gnu89 -Wextra

bmp.o: bmp.c bmp.h
	gcc -c bmp.c -obmp.o -lm

clean:
	rm bmp.o
	rm ip_lib.o
	rm main_iplib.o
	rm main_iplib
