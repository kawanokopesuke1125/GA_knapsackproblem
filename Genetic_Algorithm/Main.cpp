/*******************************************************************************
	Main.cpp
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Main.h"
#include "GA.h"


/// メイン関数
int main( int argc, char **argv )
{
	/// GA初期設定
	SetupGA();

	/// 評価値の初期化
	InitFitness();

	/// シミュレーション
	Simulation();

	return ( 0 );
}
/*******************************************************************************
	END
*******************************************************************************/
