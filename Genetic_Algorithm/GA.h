/*******************************************************************************
	GA.h
*******************************************************************************/
#ifndef ___GA___
#define ___GA___

/// GAの設定
#define NUM_OF_OBJECT (20)
#define MAX_WEIGHT (150)

#define MAX_POP_SIZE	(8)								/// 1世代あたりの遺伝子数（2以上の偶数で指定）
#define MAX_GENERATION	(100)					/// シミュレーション世代
#define NUM_OF_PARAM	(NUM_OF_OBJECT)		/// パラメータ数（最適化するパラメータ数）
#define BIT	(1)													/// 1パラメータあたりのビット数
#define TOTAL_BIT		(NUM_OF_PARAM*BIT)	/// 1個の遺伝子の総ビット数
#define POP_SIZE		(MAX_POP_SIZE)				/// 遺伝子数
#define ELITE_SIZE		(2)									/// エリート選択する遺伝子数
#define CROSS_RATE		(0.8)							/// 交差確率
#define MUTATE_RATE		(0.3)						/// 突然変異確率
#define NUM_OF_FITNESS	(1)							/// 評価項目数


/// 変数宣言
typedef struct {
	unsigned int genotype[ TOTAL_BIT +2 ];	/// 遺伝子型
	double final_fitness;										/// 評価値
} INDIVIDUAL;

typedef struct {
	unsigned int value;		/// 品物の価値
	unsigned int weight;	/// 品物の重さ
} OBJECT;


/// 関数宣言
void  SetupGA( void );
void InitFitness( void );
void GeneticOperation( void );
void WriteGtype( const unsigned int &simgene );
void SaveAllData( const unsigned int &simgene );
void SaveBestData( const unsigned int &simgene );
void Simulation( void );

#endif
/*******************************************************************************
	END
*******************************************************************************/
