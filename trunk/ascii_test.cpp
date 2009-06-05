#include <stdio.h>

int main(int argc, char *argv[])
{
  for (unsigned char c = 0; c < 256; ++c) {
    printf("%d [%c]\t", c, c);
  }
}
