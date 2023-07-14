CC = nvcc
CFLAGS = -O3 -Xcompiler -Wall
LDFLAGS = -lpng -lm

hue_modify: hue_modify.cu
	$(CC) $(CFLAGS) -o hue_modify hue_modify.cu $(LDFLAGS)

.PHONY: clean
clean:
	rm -f hue_modify
