/*******************************************************************************
	GA.cpp
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "Main.h"
#include "GA.h"

/// �O���[�o���ϐ�
static unsigned int total_simgene = 0;
static INDIVIDUAL indiv_org[ MAX_POP_SIZE +2 ];
static INDIVIDUAL indiv_cpy[ MAX_POP_SIZE +2 ];
static OBJECT object[NUM_OF_PARAM +2];


/// �O����P�̃����_���Ȑ����idouble�F�����_�ȉ�6���܂Łj
static double rand_num_double( void )
{
	int itmp = 0;
	double dtmp = 0.0;

	dtmp = rand()/((double)RAND_MAX+1.0);
	itmp = (int)(dtmp*1000000.0);
	dtmp = (double)itmp / 1000000.0;

	return ( dtmp );
}


/// �O������������̃����_���Ȑ����iunsigned int�j
static unsigned int rand_num_int( const unsigned int &max_num )
{
	return ( (unsigned int)( ( rand()/( (double)RAND_MAX+1.0 ) )*max_num ) );
}


/// �x�X�g�f�[�^�̕ۑ�
void SaveBestData( const unsigned int &simgene )
{
	static FILE *fp;

	/// �t�@�C�����I�[�v��
	if( ( fopen_s( &fp, "./data/best_data.txt", "a" ) )!=NULL ) {
		perror( NULL );
		exit( 1 );
	}
	else {
		/// ����
		fprintf(fp, "%u\t", simgene);

		/// �K���x�i�x�X�g�j
		unsigned int best_indiv = 0;
		double best_fitness = 0.0;

		best_fitness = indiv_org[0].final_fitness;
		for( int i = 0; i < MAX_POP_SIZE; i++ ) {
			if( best_fitness < indiv_org[i].final_fitness ) {
				best_indiv = i;
				best_fitness = indiv_org[i].final_fitness;
			}
			else {}
		}
		fprintf(fp, "%u\t", best_indiv);
		fprintf(fp, "%lf\t", best_fitness);

		/// �K���x�i���ρj
		double sum_fitness = 0.0;
		for( int i = 0; i < MAX_POP_SIZE; i++ ) {
			sum_fitness += indiv_org[i].final_fitness;
		}
		sum_fitness /= (double)MAX_POP_SIZE;
		fprintf(fp, "%lf\t", sum_fitness);

		/// ��`�q�^
		for( int i = 0; i < TOTAL_BIT; i++ ) {
			fprintf(fp, "%u", indiv_org[best_indiv].genotype[i]);
		}
		fprintf(fp, "\n");

		/// �t�@�C�����N���[�Y
		fclose( fp );
	}
}


/// �S��`�q�f�[�^�̕ۑ�
void SaveAllData( const unsigned int &simgene )
{
	static FILE *fp;

	/// �t�@�C�����I�[�v��
	if( ( fopen_s( &fp, "./data/all_data.txt", "a" ) )!=NULL ) {
		perror( NULL );
		exit( 1 );
	}
	else {
		for( int i = 0; i < MAX_POP_SIZE; i++ ) {
			/// ����
			fprintf(fp, "%u\t", simgene);

			/// ��`�q�ԍ�
			fprintf(fp, "%d\t", i);

			/// �]���l
			fprintf(fp, "%lf\t", indiv_org[i].final_fitness);

			/// ��`�q�^
			for( int j = 0; j < TOTAL_BIT; j++ ) {
				fprintf(fp, "%u", indiv_org[i].genotype[j]);
			}
			fprintf(fp, "\n");
		}

		/// �t�@�C�����N���[�Y
		fclose( fp );
	}
}


/// �Z�[�u�t�@�C���̍쐬
static void MakeFileIndex( void )
{
	static FILE *fp;

	/// �t�@�C�����I�[�v���iall_data.txt�j
	if( (fopen_s( &fp, "./data/all_data.txt", "a" ))!=NULL ) {
		perror( NULL );
		exit( 1 );
	}
	else {
		struct tm newtime;
		__time64_t long_time;
		errno_t err;

		_time64( &long_time );
		err = _localtime64_s( &newtime, &long_time);
		if( err ) {
			perror( NULL );
			exit( 1 );
		}
		else {
			fprintf(fp, "\n");
			fprintf(fp, "%2d/%2d/%2d, %2d:%2d:%2d \n", newtime.tm_year+1900, newtime.tm_mon+1, newtime.tm_mday, newtime.tm_hour, newtime.tm_min, newtime.tm_sec);
		}

		/// �C���f�b�N�X
		/// ����
		fprintf(fp, "Generation\t");

		/// ��`�q�ԍ�
		fprintf(fp, "Indiv No.\t");

		/// �]���l
		fprintf(fp, "Fitness\t");

		/// ��`�q�^
		fprintf(fp, "Genotype\t");

		fprintf(fp, "\n");

		/// �t�@�C�����N���[�Y
		fclose( fp );
	}

	/// �t�@�C�����I�[�v���ibest_data.txt�j
	if( ( fopen_s( &fp, "./data/best_data.txt", "a" ) )!=NULL ) {
		perror( NULL );
		exit( 1 );
	}
	else {
		struct tm newtime;
		__time64_t long_time;
		errno_t err;

		_time64( &long_time );
		err = _localtime64_s( &newtime, &long_time);
		if( err ) {
			perror( NULL );
			exit( 1 );
		}
		else {
			fprintf(fp, "\n");
			fprintf(fp, "%2d/%2d/%2d, %2d:%2d:%2d \n", newtime.tm_year+1900, newtime.tm_mon+1, newtime.tm_mday, newtime.tm_hour, newtime.tm_min, newtime.tm_sec);
		}

		/// �C���f�b�N�X
		/// ����
		fprintf(fp, "Generation\t");

		/// ��`�q�ԍ�
		fprintf(fp, "Indiv No.\t");

		/// �]���l
		fprintf(fp, "Fitness(Best)\t");
		fprintf(fp, "Fitness(Ave.)\t");

		/// ��`�q�^
		fprintf(fp, "Genotype\t");

		fprintf(fp, "\n");

		/// �t�@�C�����N���[�Y
		fclose( fp );
	}
}


/// ���ꂩ��V�~�����[�V���������`�q����ۑ�
void WriteGtype( const unsigned int &simgene )
{
	static FILE *fp;

	/// �t�@�C�����I�[�v��
	if( ( fopen_s( &fp, "./data/gtype.txt", "w" ) )!=NULL ) {
		perror( NULL );
		exit( 1 );
	}
	else {
		/// ���V�~�����[�V�������㐔�̏����o��
		fprintf(fp, "%u \n", total_simgene+simgene+1);

		/// �f�[�^�̏����o��
		for( int i = 0; i < MAX_POP_SIZE; i++ ) {
			for( int j = 0; j < TOTAL_BIT; j++ ) {
				fprintf( fp, "%u ", indiv_org[i].genotype[j] );
			}
			fprintf( fp, "\n" );
		}
		fprintf( fp, "\n" );

		/// �t�@�C�����N���[�Y
		fclose( fp );
	}
}


/// �ˑR�ψ�
static void Mutation( void )
{
	/// �ˑR�ψق́u1�v����X�^�[�g�i�[���Ԗڂ̃G���[�g�͂��̂܂ܕۑ��j
	for( int i = 1; i < MAX_POP_SIZE; i++ ) {
		//printf("\t mutate : %3d \t", i);
		for( int j = 0; j < TOTAL_BIT; j++ ) {
			if( rand_num_double() < MUTATE_RATE ) {
				if( indiv_cpy[i].genotype[j] != 0 ) {
					indiv_cpy[i].genotype[j] = 0;
				}
				else {
					indiv_cpy[i].genotype[j] = 1;
				}
				//printf("%d, ", j);
			}
			else {}
		}
		//printf("\n");
	}
}

/// ���[���b�g�I��
static unsigned int Roulette( void )
{
	/// ��`�q�Q�̕]���l�̑��a���v�Z
	double sum_fitness = 0.0;
	for( int i = 0; i < MAX_POP_SIZE; i++ ) {
		sum_fitness += indiv_org[i].final_fitness*indiv_org[i].final_fitness;
	}

	/// 臒l������
	double threshold = sum_fitness * rand_num_double();

	/// ���[���b�g�I��
	unsigned int j = 0;
	double part_sum_fitness = 0.0;

	while( 1 ) {
		if( !(j < MAX_POP_SIZE) ) {
			j = MAX_POP_SIZE-1;
			break;
		}
		else {
			part_sum_fitness += indiv_org[j].final_fitness*indiv_org[j].final_fitness;
			if( part_sum_fitness < threshold ) {
				j++;
			}
			else {
				break;
			}
		}
	}

	return ( j );
}


/// 2�_����
static void Crossover( void )
{
	for( int i = ELITE_SIZE; i < MAX_POP_SIZE; i+=2 ) {
		/// ��`�q�̑I��
		unsigned int iID1 = Roulette();
		unsigned int iID2 = Roulette();

		/// �Ȃ�ׂ��قȂ��`�q���I������邽�߂̏���
		while( iID1 == iID2 ) {
			static unsigned int itmp = 0;
			if( itmp > MAX_POP_SIZE/2 ) {
				printf("\t warning: same populations are selected. \n");
				break;
			}
			else {
				itmp++;
			}
			iID2 = Roulette();
		}

		/// �����|�C���g�̌���
		unsigned int locus1 = rand_num_int( TOTAL_BIT );
		unsigned int locus2 = rand_num_int( TOTAL_BIT );

		/// �����|�C���g��locus1 < locus2�Ƃ���
		if( locus1 > locus2 ) {
			unsigned int itmp = locus2;
			locus2 = locus1;
			locus1 = itmp;
		}
		else {}

		/// 2�_�������ĕۑ�
		if( rand_num_double() < CROSS_RATE ) {
			for( int j = 0; j < (signed)locus1; j++ ) {
				indiv_cpy[i  ].genotype[j] = indiv_org[iID1].genotype[j];
				indiv_cpy[i+1].genotype[j] = indiv_org[iID2].genotype[j];
			}
			for( int j = locus1; j < (signed)locus2; j++ ) {
				indiv_cpy[i  ].genotype[j] = indiv_org[iID2].genotype[j];
				indiv_cpy[i+1].genotype[j] = indiv_org[iID1].genotype[j];
			}
			for( int j = locus2; j < TOTAL_BIT; j++ ) {
				indiv_cpy[i  ].genotype[j] = indiv_org[iID1].genotype[j];
				indiv_cpy[i+1].genotype[j] = indiv_org[iID2].genotype[j];
			}
			indiv_cpy[i  ].genotype[TOTAL_BIT] = '\0';
			indiv_cpy[i+1].genotype[TOTAL_BIT] = '\0';

			printf("\t cross  : %3u, %3u \t %3u %3u \n", iID1, iID2, locus1, locus2);
		}
		else {
			for( int j = 0; j < TOTAL_BIT; j++ ) {
				indiv_cpy[i  ].genotype[j] = indiv_org[iID1].genotype[j];
				indiv_cpy[i+1].genotype[j] = indiv_org[iID2].genotype[j];
			}
			indiv_cpy[i  ].genotype[TOTAL_BIT] = '\0';
			indiv_cpy[i+1].genotype[TOTAL_BIT] = '\0';
			printf("\t cross  : %3u, %3u \t no cross \n", iID1, iID2);
		}
	}
}


/// �G���[�g�I��
static void EliteSelection( void )
{
	unsigned int best_indiv = 0;
	double best_fitness = 0.0;

	/// �G���[�g��`�q�̒T��
	best_fitness = indiv_org[0].final_fitness;
	for( int i = 0; i < MAX_POP_SIZE; i++ ) {
		if( best_fitness < indiv_org[i].final_fitness ) {
			best_indiv = i;
			best_fitness = indiv_org[i].final_fitness;
		}
		else {}
	}
	printf("\t elite  : %3u \n", best_indiv);

	/// �G���[�g��`�q�̕ۑ�
	for( int i = 0; i < ELITE_SIZE; i++ ) {
		for( int j = 0; j < TOTAL_BIT; j++ ) {
			indiv_cpy[i].genotype[j] = indiv_org[best_indiv].genotype[j];
		}
	}
}


/// ��`�q����
void GeneticOperation( void )
{
	EliteSelection();
	Crossover();
	Mutation();

	/// ��`�q�̃R�s�[��������
	for( int i = 0; i < MAX_POP_SIZE; i++ ) {
		for( int j = 0; j < TOTAL_BIT; j++ ) {
			indiv_org[i].genotype[j] = indiv_cpy[i].genotype[j];
		}
		indiv_org[i].final_fitness = 0.0;
	}
}


/// �]���l�̏�����
void InitFitness( void )
{
	/// �]���l�̏�����
	for( int i = 0; i < MAX_POP_SIZE+2; i++ ) {
		indiv_org[i].final_fitness = 0.0;
		indiv_cpy[i].final_fitness = 0.0;
	}
}


/// GA�����ݒ�
void SetupGA( void )
{
	/// �����V�[�h�̏�����
	srand( (unsigned int)time(NULL) );

	/// ������`�q�̐���
	for( int i = 0; i < MAX_POP_SIZE; i++ ) {
		for( int j = 0; j < TOTAL_BIT; j++ ) {
			indiv_org[i].genotype[j] = rand_num_int( 2 );
			indiv_cpy[i].genotype[j] = indiv_org[i].genotype[j];
		}
	}

	/// �Z�[�u�t�@�C���̍쐬
	MakeFileIndex();

	/// �i�b�v�T�b�N�̐ݒ�
	for( int i = 0; i < NUM_OF_OBJECT; i++ ) {
		object[i].value = 0;
		object[i].weight = 0;
	}
	object[0].value   = 21;
	object[0].weight = 2;
	object[1].value   = 22;
	object[1].weight = 10;
	object[2].value   = 28;
	object[2].weight = 7;
	object[3].value   = 21;
	object[3].weight = 2;
	object[4].value   = 12;
	object[4].weight = 4;
	object[5].value   = 24;
	object[5].weight = 9;
	object[6].value   = 15;
	object[6].weight = 10;
	object[7].value   = 2;
	object[7].weight = 7;
	object[8].value   = 25;
	object[8].weight = 8;
	object[9].value   = 28;
	object[9].weight = 5;

	object[10].value   = 4;
	object[10].weight = 3;
	object[11].value   = 22;
	object[11].weight = 10;
	object[12].value   = 36;
	object[12].weight = 9;
	object[13].value   = 2;
	object[13].weight = 8;
	object[14].value   = 7;
	object[14].weight = 8;
	object[15].value   = 40;
	object[15].weight = 5;
	object[16].value   = 14;
	object[16].weight = 7;
	object[17].value   = 40;
	object[17].weight = 3;
	object[18].value   = 33;
	object[18].weight = 9;
	object[19].value   = 21;
	object[19].weight = 7;

	object[20].value   = 28;
	object[20].weight = 2;
	/*/object[21].value   = 22;
	object[21].weight = 10;
	object[22].value   = 14;
	object[22].weight = 7;
	object[23].value   = 36;
	object[23].weight = 9;
	object[24].value   = 28;
	object[24].weight = 7;
	object[25].value   = 21;
	object[25].weight = 2;
	object[26].value   = 18;
	object[26].weight = 10;
	object[27].value   = 12;
	object[27].weight = 4;
	object[28].value   = 24;
	object[28].weight = 9;
	object[29].value   = 15;
	object[29].weight = 10;

	object[30].value   = 21;
	object[30].weight = 4;
	object[31].value   = 2;
	object[31].weight = 7;
	object[32].value   = 25;
	object[32].weight = 8;
	object[33].value   = 28;
	object[33].weight = 5;
	object[34].value   = 28;
	object[34].weight = 2;
	object[35].value   = 4;
	object[35].weight = 3;
	object[36].value   = 22;
	object[36].weight = 10;
	object[37].value   = 36;
	object[37].weight = 9;
	object[38].value   = 31;
	object[38].weight = 7;
	object[39].value   = 2;
	object[39].weight = 8;

	object[40].value   = 7;
	object[40].weight = 8;
	object[41].value   = 40;
	object[41].weight = 5;
	object[42].value   = 14;
	object[42].weight = 7;
	object[43].value   = 4;
	object[43].weight = 5;
	object[44].value   = 28;
	object[44].weight = 7;
	object[45].value   = 40;
	object[45].weight = 3;
	object[46].value   = 33;
	object[46].weight = 9;
	object[47].value   = 35;
	object[47].weight = 7;
	object[48].value   = 21;
	object[48].weight = 7;
	object[49].value   = 20;
	object[49].weight = 9;
	/*/
	/*
	for( int i = 0; i < NUM_OF_OBJECT; i++ ) {
		printf("%d\t", object[i].value);
		printf("%d\n", object[i].weight);
	}
	system("pause");
	*/

	/*
	int itmp_sum = 0;
	for( int i = 0; i < NUM_OF_OBJECT; i++ ) {
		itmp_sum += object[i].weight;
	}
	printf("total weight = %d\n", itmp_sum);
	system("pause");
	*/
}


/// �]��
void Simulation( void )
{
	static unsigned int simindiv = 0;		/// �V�~�����[�V������
	static unsigned int simgene = 0;		/// �V�~�����[�V��������

	static unsigned int calc_weight = 0;

	while(1) {
		/// �o�߂�\��
		printf("simgene: %3u, simindiv: %3u, \n", simgene, simindiv);

		/// 1�̌́isimindiv�Ԗځj��]������
		calc_weight = 0;
		for( int j = 0; j < NUM_OF_OBJECT; j++ ) {
			calc_weight += indiv_org[simindiv].genotype[j] * object[j].weight;
		}
		//printf("calc_weight = %d\n", calc_weight);	
		//system("pause");

		if( calc_weight > MAX_WEIGHT ) {
			indiv_org[simindiv].final_fitness = 0.0;
		}
		else {
			for( int j = 0; j < NUM_OF_OBJECT; j++ ) {
				indiv_org[simindiv].final_fitness += (double)indiv_org[simindiv].genotype[j] * object[j].value;
			}
		}

		simindiv++;		/// ���̌̂ɃJ�E���g�A�b�v
		if( simindiv > MAX_POP_SIZE-1 ) {
			/// 1����̕]�����I�������
			simindiv = 0;
			/// ��`�q�f�[�^�̕ۑ�
			SaveAllData( simgene );
			SaveBestData( simgene );

			simgene++;		/// ���̐���ɃJ�E���g�A�b�v
			if( simgene > MAX_GENERATION-1 ) {
				/// �S����̃V�~�����[�V�����I��
				break;
			}
			else {
				/// ��`�q����
				GeneticOperation();

				/// �]���l�̏������i�V�����������`�q�̕]���l���[���ɂ���j
				InitFitness();

				/// ���ꂩ��V�~�����[�V���������`�q����ۑ�
				WriteGtype( simgene );
			}
		}
		else {}
	}
}
/*******************************************************************************
	END
*******************************************************************************/
