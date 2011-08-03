#include <stdio.h>
#include <stdlib.h>
#include "data.h"



int IsShoreline(int y, int x, TOPOPIX ** TopoMap){
	if(TopoMap[y][x].Shoreline){
			return 1;
	}
	else return 0;
}