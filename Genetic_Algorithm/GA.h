/*******************************************************************************
	GA.h
*******************************************************************************/
#ifndef ___GA___
#define ___GA___

/// GA�̐ݒ�
#define NUM_OF_OBJECT (20)
#define MAX_WEIGHT (150)

#define MAX_POP_SIZE	(8)								/// 1���゠����̈�`�q���i2�ȏ�̋����Ŏw��j
#define MAX_GENERATION	(100)					/// �V�~�����[�V��������
#define NUM_OF_PARAM	(NUM_OF_OBJECT)		/// �p�����[�^���i�œK������p�����[�^���j
#define BIT	(1)													/// 1�p�����[�^������̃r�b�g��
#define TOTAL_BIT		(NUM_OF_PARAM*BIT)	/// 1�̈�`�q�̑��r�b�g��
#define POP_SIZE		(MAX_POP_SIZE)				/// ��`�q��
#define ELITE_SIZE		(2)									/// �G���[�g�I�������`�q��
#define CROSS_RATE		(0.8)							/// �����m��
#define MUTATE_RATE		(0.3)						/// �ˑR�ψيm��
#define NUM_OF_FITNESS	(1)							/// �]�����ڐ�


/// �ϐ��錾
typedef struct {
	unsigned int genotype[ TOTAL_BIT +2 ];	/// ��`�q�^
	double final_fitness;										/// �]���l
} INDIVIDUAL;

typedef struct {
	unsigned int value;		/// �i���̉��l
	unsigned int weight;	/// �i���̏d��
} OBJECT;


/// �֐��錾
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
