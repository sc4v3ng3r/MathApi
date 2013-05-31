#ifndef DEFINES_H
#define DEFINES_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <errno.h>

#define TRUE 1
#define FALSE 0

#ifdef WIN32
#define CLEAR_SCREEN system("cls");
#else
#define CLEAR_SCREEN printf("\33[H\33[2J");
#endif

typedef unsigned int uint;
typedef unsigned long ulong;
typedef unsigned short ushort;
typedef ushort BOOL;
#endif