// This code requires two files: keys.txt and startingpoint.txt
// Random 100 keys are stored in "keys.txt".
// Number of data of an experiment can be increased after an experiment is completed by storing the bias results in the "startingpoint.txt" file.
#include <Windows.h>
#include <stdio.h>
#include <cuda.h>
#include <time.h>
typedef unsigned __int64 bit64;
typedef unsigned long bit32;
typedef unsigned char bit8;
#define ROL(x,r) ((x) = ((x) << (r)) | ((x) >> (32-(r))))
double PCFreq = 0.0;
__int64 CounterStart = 0;
__int64 blocksize = 512, threadsize = 1024, loop = 1024 * 256;
__int64 totalthreads = blocksize * threadsize;
//#define loop_d 128*1024
#define loop_d 1024*128
bit32 plaintext[4], ciphertext[4], key[8], w[132], k[132], phi = 0x9e3779b9;
bit32 x[4]; //plaintext words
int difference, devicenumber,userloop,pairs=1,numberofrounds, startingpoint=0;
float experiment=1; // This value should be an integer in [1,100] and determines the number of repeatition of the experiment with different keys
// key.txt contains 100 random keys (as 400 lines containing random 32-bit values in each line)

/* Serpent's S-boxes */
bit32 S[8][16] = {
	{ 0x3, 0x8, 0xf, 0x1, 0xa, 0x6, 0x5, 0xb, 0xe, 0xd, 0x4, 0x2, 0x7, 0x0, 0x9, 0xc },
	{ 0xf, 0xc, 0x2, 0x7, 0x9, 0x0, 0x5, 0xa, 0x1, 0xb, 0xe, 0x8, 0x6, 0xd, 0x3, 0x4 },
	{ 0x8, 0x6, 0x7, 0x9, 0x3, 0xc, 0xa, 0xf, 0xd, 0x1, 0xe, 0x4, 0x0, 0xb, 0x5, 0x2 },
	{ 0x0, 0xf, 0xb, 0x8, 0xc, 0x9, 0x6, 0x3, 0xd, 0x1, 0x2, 0x4, 0xa, 0x7, 0x5, 0xe },
	{ 0x1, 0xf, 0x8, 0x3, 0xc, 0x0, 0xb, 0x6, 0x2, 0x5, 0x4, 0xa, 0x9, 0xe, 0x7, 0xd },
	{ 0xf, 0x5, 0x2, 0xb, 0x4, 0xa, 0x9, 0xc, 0x0, 0x3, 0xe, 0x8, 0xd, 0x6, 0x7, 0x1 },
	{ 0x7, 0x2, 0xc, 0x5, 0x8, 0x4, 0x6, 0xb, 0xe, 0x9, 0x1, 0xf, 0xd, 0x3, 0xa, 0x0 },
	{ 0x1, 0xd, 0xf, 0x0, 0xe, 0x8, 0x2, 0xb, 0x7, 0x4, 0xc, 0xa, 0x9, 0x3, 0x5, 0x6 }
};
void StartCounter(){
	LARGE_INTEGER li;
	if (!QueryPerformanceFrequency(&li))
		printf("QueryPerformanceFrequency failed!\n");

	PCFreq = double(li.QuadPart) / 1000.0;

	QueryPerformanceCounter(&li);
	CounterStart = li.QuadPart;
}
double GetCounter(){
	LARGE_INTEGER li;
	QueryPerformanceCounter(&li);
	return double(li.QuadPart - CounterStart) / PCFreq;
}
void get_user_inputs() {
	key[4] = 0x80000000; key[5] = 0x00000000; key[6] = 0x00000000; key[7] = 0x00000000;
	printf("Select Cuda Device (0,1,2): ");
	scanf("%d", &devicenumber);
	printf("Choose an experiment:\n");
	printf(
		"(3) 3-round DL Experiment for [PZWD24]\n"
		"(32) SERPENT Benchmark\n"
		"(4) 4-round DL Experiment for our distinguisher (14 bit left rotation of [DIK08])\n"
		"(5) 5-round DL Experiment for [HDE24]\n"
		"(66) 6-round DL Experiment for [PZWD24] (final correction at eprint)\n"
		"(74) 4-round DL Experiment: Middle 4 rounds of 7-round DL of [PZWD24] \n"
		"(75) 5-round DL Experiment: Middle 5 rounds of 7-round DL of [PZWD24] \n"
		"(76) 6-round DL Experiment: Last 6 rounds of 7-round DL of [PZWD24] \n"
		"Choice: "
	);
	scanf("%d", &numberofrounds);
	printf("Number of Pairs: 2^36+");
	scanf("%d", &userloop);
	printf("Starting point (default 0): ");
	scanf("%d", &startingpoint);
	for (int i = 0; i<userloop; i++) pairs *= 2;
}
void key_s(int counter, int sbox) {
	int i, j;
	bit32 temp, x_temp[4];
	for (i = 0; i<4; i++) x_temp[i] = 0;
	for (i = 0; i<32; i++) {
		temp = 0;
		for (j = 0; j<4; j++) temp = temp | (((w[4 * counter + j] >> i) & 0x1) << j);
		temp = S[sbox][temp];
		for (j = 0; j<4; j++) x_temp[j] = x_temp[j] | (((temp >> j) & 0x1) << i);
	}
	for (i = 0; i<4; i++) k[4 * counter + i] = x_temp[i];
}
void key_schedule() {
	int i;
	for (i = 0; i<8; i++) w[i] = key[i];
	for (i = 8; i<16; i++) {
		w[i] = w[i - 8] ^ w[i - 5] ^ w[i - 3] ^ w[i - 1] ^ phi ^ (i - 8);
		w[i] = (w[i] << 11) | (w[i] >> 21);
	}
	for (i = 0; i<8; i++) w[i] = w[i + 8];
	// Generate w[i]'s
	for (i = 8; i<132; i++) {
		w[i] = (w[i - 8] ^ w[i - 5] ^ w[i - 3] ^ w[i - 1] ^ phi^i);
		// w[i] = w[i] <<< 11
		w[i] = (w[i] << 11) | (w[i] >> 21);
	}
	// Generate k[i]'s
	key_s(0, 3);
	key_s(1, 2);	key_s(2, 1);	key_s(3, 0);	key_s(4, 7);	key_s(5, 6);	key_s(6, 5);	key_s(7, 4);
	key_s(8, 3);	key_s(9, 2);	key_s(10, 1); key_s(11, 0); key_s(12, 7); key_s(13, 6); key_s(14, 5); key_s(15, 4);
	key_s(16, 3); key_s(17, 2); key_s(18, 1); key_s(19, 0); key_s(20, 7); key_s(21, 6); key_s(22, 5); key_s(23, 4);
	key_s(24, 3); key_s(25, 2); key_s(26, 1); key_s(27, 0); key_s(28, 7); key_s(29, 6); key_s(30, 5); key_s(31, 4); key_s(32, 3);
}
void s(int sbox) {
	int i, j;
	bit32 temp, x_temp[4];
	for (i = 0; i<4; i++) x_temp[i] = 0;
	for (i = 0; i<32; i++) {
		temp = 0;
		for (j = 0; j<4; j++) temp = temp | (((x[j] >> i) & 0x1) << j);
		temp = S[sbox][temp];
		for (j = 0; j<4; j++) x_temp[j] = x_temp[j] | (((temp >> j) & 0x1) << i);
	}
	for (i = 0; i<4; i++) x[i] = x_temp[i];
}
void Sb0(bit32 a, bit32 b, bit32 c, bit32 d) {
	bit32 t1, t3, t4, t7, t12;
	t1 = a ^ d;
	t3 = c ^ t1;
	t4 = b ^ t3;
	x[3] = (a & d) ^ t4;
	t7 = a ^ (b & t1);
	x[2] = t4 ^ (c | t7);
	t12 = x[3] & (t3 ^ t7);
	x[1] = (~t3) ^ t12;
	x[0] = t12 ^ (~t7);
}
void Sb1(bit32 a, bit32 b, bit32 c, bit32 d) {
	bit32 t2, t5, t7, t8, t11;
	t2 = b ^ (~a);
	t5 = c ^ (a | t2);
	x[2] = d ^ t5;
	t7 = b ^ (d | t2);
	t8 = t2 ^ x[2];
	x[3] = t8 ^ (t5 & t7);
	t11 = t5 ^ t7;
	x[1] = x[3] ^ t11;
	x[0] = t5 ^ (t8 & t11);
}
void Sb2(bit32 a, bit32 b, bit32 c, bit32 d)   {
	bit32 t1, t2, t3, t5, t6, t7;
	t1 = ~a;
	t2 = b ^ d;
	t3 = c & t1;
	x[0] = t2 ^ t3;
	t5 = c ^ t1;
	t6 = c ^ x[0];
	t7 = b & t6;
	x[3] = t5 ^ t7;
	x[2] = a ^ ((d | t7) & (x[0] | t5));
	x[1] = (t2 ^ x[3]) ^ (x[2] ^ (d | t1));
}
void Sb3(bit32 a, bit32 b, bit32 c, bit32 d) {
	bit32 t1, t2, t3, t4, t5, t6, t8, t9, t10, t12;
	t1 = a ^ b;
	t2 = a & c;
	t3 = a | d;
	t4 = c ^ d;
	t5 = t1 & t3;
	t6 = t2 | t5;
	x[2] = t4 ^ t6;
	t8 = b ^ t3;
	t9 = t6 ^ t8;
	t10 = t4 & t9;
	x[0] = t1 ^ t10;
	t12 = x[2] & x[0];
	x[1] = t9 ^ t12;
	x[3] = (b | d) ^ (t4 ^ t12);
}
void Sb4(bit32 a, bit32 b, bit32 c, bit32 d){
	bit32 t1, t2, t3, t4, t6, t7, t9, t10, t11;
	t1 = a ^ d;
	t2 = d & t1;
	t3 = c ^ t2;
	t4 = b | t3;
	x[3] = t1 ^ t4;
	t6 = ~b;
	t7 = t1 | t6;
	x[0] = t3 ^ t7;
	t9 = a & x[0];
	t10 = t1 ^ t6;
	t11 = t4 & t10;
	x[2] = t9 ^ t11;
	x[1] = (a ^ t3) ^ (t10 & x[2]);
}
void Sb5(bit32 a, bit32 b, bit32 c, bit32 d) {
	bit32 t1, t2, t3, t4, t5, t7, t8, t10, t11, t12;
	t1 = ~a;
	t2 = a ^ b;
	t3 = a ^ d;
	t4 = c ^ t1;
	t5 = t2 | t3;
	x[0] = t4 ^ t5;
	t7 = d & x[0];
	t8 = t2 ^ x[0];
	x[1] = t7 ^ t8;
	t10 = t1 | x[0];
	t11 = t2 | t7;
	t12 = t3 ^ t10;
	x[2] = t11 ^ t12;
	x[3] = (b ^ t7) ^ (x[1] & t12);
}
void Sb6(bit32 a, bit32 b, bit32 c, bit32 d) {
	bit32 t1, t2, t3, t4, t5, t7, t8, t9, t11;
	t1 = ~a;
	t2 = a ^ d;
	t3 = b ^ t2;
	t4 = t1 | t2;
	t5 = c ^ t4;
	x[1] = b ^ t5;
	t7 = t2 | x[1];
	t8 = d ^ t7;
	t9 = t5 & t8;
	x[2] = t3 ^ t9;
	t11 = t5 ^ t8;
	x[0] = x[2] ^ t11;
	x[3] = (~t5) ^ (t3 & t11);
}
void Sb7(bit32 a, bit32 b, bit32 c, bit32 d) {
	bit32 t1, t2, t3, t4, t5, t6, t8, t9, t11, t12;
	t1 = b ^ c;
	t2 = c & t1;
	t3 = d ^ t2;
	t4 = a ^ t3;
	t5 = d | t1;
	t6 = t4 & t5;
	x[1] = b ^ t6;
	t8 = t3 | x[1];
	t9 = a & t4;
	x[3] = t1 ^ t9;
	t11 = t4 ^ t8;
	t12 = x[3] & t11;
	x[2] = t3 ^ t12;
	x[0] = (~t11) ^ (x[3] & x[2]);
}
void key_addition(int round) {	for (int i = 0; i<4; i++) x[i] = x[i] ^ k[4 * round + i];}
void linear_transformation() {
	// x0 = x0 <<< 13
	x[0] = (x[0] << 13) | x[0] >> 19;
	// x2 <<< 3
	x[2] = (x[2] << 3) | (x[2] >> 29);
	// x1 = x1^x0^x2
	x[1] = x[1] ^ x[0] ^ x[2];
	// x3 = x3^x2^(x0 << 3)
	x[3] = x[3] ^ x[2] ^ (x[0] << 3);
	// x1 = x1 <<< 1
	x[1] = (x[1] << 1) | (x[1] >> 31);
	// x3 = x3 <<< 7
	x[3] = (x[3] << 7) | (x[3] >> 25);
	// x0 = x0 ^ x1 ^ x3
	x[0] = x[0] ^ x[1] ^ x[3];
	// x2 = x2 ^ x3 ^ (x1 << 7)
	x[2] = x[2] ^ x[3] ^ (x[1] << 7);
	// x0 = x0 <<< 5
	x[0] = (x[0] << 5) | (x[0] >> 27);
	// x2 = x2 <<< 22
	x[2] = (x[2] << 22) | (x[2] >> 10);
}
__global__ void Serpent4round(bit32 p0, bit32 p1, bit32 p2, bit32 p3, bit64 *hit, bit32 *k_d) {
	__shared__ bit32  k[132];
	bit32 x0, x1, x2, x3, plaintext0 = p0, plaintext1 = p1, plaintext2 = p2, plaintext3 = p3 + blockIdx.x *blockDim.x + threadIdx.x;
	bit32 a, b, c, d;
	bit32 ciphertext0, ciphertext1, ciphertext2, ciphertext3;
	bit32 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	bit32 round;
	bit64 counter = 0;
	if (threadIdx.x<132) k[threadIdx.x] = k_d[threadIdx.x];
	__syncthreads();
	for (int i = 0; i < loop_d; i++) {
		x0 = plaintext0; x1 = plaintext1; x2 = plaintext2; x3 = plaintext3;
		// Round 2
		round = 2;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 3
		round = 3;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 4
		round = 4;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		round = 5;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		ciphertext0 = x0; ciphertext1 = x1; ciphertext2 = x2; ciphertext3 = x3;

		// Give difference
		x0 = plaintext0 ^ (0x00000010<<14);
		x1 = plaintext1;
		x2 = plaintext2 ^ (0x00000090<<14);
		x3 = plaintext3;
		// Round 2
		round = 2;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 3
		round = 3;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 4
		round = 4;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		round = 5;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);

		// Check Difference
		x2 ^= ciphertext2;
		x3 ^= ciphertext3;
		t1 = (x2 >> 29) & 0x1;
		t2 = x3 & 0x1;
		t1 ^= t2;
		if (t1 == 0) counter++;

		// Change plaintext
		plaintext0 = x0; plaintext1 = x1; plaintext2 = x2; plaintext3 = x3;
	}
	hit[blockIdx.x *blockDim.x + threadIdx.x] = counter;
}
__global__ void Serpent3round_PZWD24(bit32 p0, bit32 p1, bit32 p2, bit32 p3, bit64* hit, bit32* k_d) {
	__shared__ bit32  k[132];
	bit32 x0, x1, x2, x3, plaintext0 = p0, plaintext1 = p1, plaintext2 = p2, plaintext3 = p3 + blockIdx.x * blockDim.x + threadIdx.x;
	bit32 a, b, c, d;
	bit32 ciphertext0, ciphertext1, ciphertext2, ciphertext3;
	bit32 t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	bit32 round;
	bit64 counter = 0;
	if (threadIdx.x < 132) k[threadIdx.x] = k_d[threadIdx.x];
	__syncthreads();
	for (int i = 0; i < loop_d; i++) {
		x0 = plaintext0; x1 = plaintext1; x2 = plaintext2; x3 = plaintext3;
		// Round 5
		round = 5;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 6
		round = 6;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 7
		round = 7;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = b ^ c; t2 = c & t1; t3 = d ^ t2; t4 = a ^ t3; t5 = d | t1; t6 = t4 & t5; x1 = b ^ t6; t8 = t3 | x1; t9 = a & t4; x3 = t1 ^ t9; t11 = t4 ^ t8; t12 = x3 & t11; x2 = t3 ^ t12; x0 = (~t11) ^ (x3 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);

		// Give difference
		ciphertext0 = x0; ciphertext1 = x1; ciphertext2 = x2; ciphertext3 = x3;		
		x0 = plaintext0 ^ 0x01000000;
		x1 = plaintext1 ^ 0x01000000;
		x2 = plaintext2 ^ 0x00000000;
		x3 = plaintext3 ^ 0x00000000;

		// Round 5
		round = 5;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 6
		round = 6;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 7
		round = 7;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = b ^ c; t2 = c & t1; t3 = d ^ t2; t4 = a ^ t3; t5 = d | t1; t6 = t4 & t5; x1 = b ^ t6; t8 = t3 | x1; t9 = a & t4; x3 = t1 ^ t9; t11 = t4 ^ t8; t12 = x3 & t11; x2 = t3 ^ t12; x0 = (~t11) ^ (x3 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);

		// Check Difference
/*		t1 = x0 ^ ciphertext0;
		t2 = x1 ^ ciphertext1;
		t3 = x2 ^ ciphertext2;
		t4 = x3 ^ ciphertext3;

		bit32 result = 0;
		t1 = t1 & 0x08010812;
		t2 = t2 & 0xa00a2004;
		t3 = t3 & 0x00000010;
		t4 = t4 & 0x000a0000;
		for (int i = 0; i < 32; i++) {
			result ^= ((t1 >> i) & 0x1);
			result ^= ((t2 >> i) & 0x1);
			result ^= ((t3 >> i) & 0x1);
			result ^= ((t4 >> i) & 0x1);
		}
		if (result == 0) counter++;*/

/*		t0 = x0 & 0x44900010;
		t1 = x1 & 0x00010000;
		t2 = x2 & 0x01980010;
		t3 = x3 & 0x12000200;*/
		x0 ^= ciphertext0;
		x1 ^= ciphertext1;
		x2 ^= ciphertext2;
		x3 ^= ciphertext3;

		t0 = x0 & 0x12000200;
		t1 = x1 & 0x01980010;
		t2 = x2 & 0x00010000;
		t3 = x3 & 0x44900010;

		t0 ^= t0 >> 1;
		t0 ^= t0 >> 2;
		t0 = (t0 & 0x11111111U) * 0x11111111U;
		t0 = (t0 >> 28) & 1;

		t1 ^= t1 >> 1;
		t1 ^= t1 >> 2;
		t1 = (t1 & 0x11111111U) * 0x11111111U;
		t1 = (t1 >> 28) & 1;

		t2 ^= t2 >> 1;
		t2 ^= t2 >> 2;
		t2 = (t2 & 0x11111111U) * 0x11111111U;
		t2 = (t2 >> 28) & 1;

		t3 ^= t3 >> 1;
		t3 ^= t3 >> 2;
		t3 = (t3 & 0x11111111U) * 0x11111111U;
		t3 = (t3 >> 28) & 1;

		t0 = t0 ^ t1 ^ t2 ^ t3;

		if (t0 == 0) counter++;

		// Change plaintext
		plaintext0 = x0; plaintext1 = x1; plaintext2 = x2; plaintext3 = x3;
	}
	hit[blockIdx.x * blockDim.x + threadIdx.x] = counter;
}
__global__ void Serpent6round_PZWD24_final(bit32 p0, bit32 p1, bit32 p2, bit32 p3, bit64* hit, bit32* k_d) {
	__shared__ bit32  k[132];
	bit32 x0, x1, x2, x3, plaintext0 = p0, plaintext1 = p1, plaintext2 = p2, plaintext3 = p3 + blockIdx.x * blockDim.x + threadIdx.x;
	bit32 a, b, c, d;
	bit32 ciphertext0, ciphertext1, ciphertext2, ciphertext3;
	bit32 t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	bit32 round;
	bit64 counter = 0;
	if (threadIdx.x < 132) k[threadIdx.x] = k_d[threadIdx.x];
	__syncthreads();
	for (int i = 0; i < loop_d; i++) {
		x0 = plaintext0; x1 = plaintext1; x2 = plaintext2; x3 = plaintext3;
		// Round 4
		x0 ^= k[16];	x1 ^= k[17];	x2 ^= k[18];	x3 ^= k[19];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		round = 5;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 6
		round = 6;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 7
		round = 7;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = b ^ c; t2 = c & t1; t3 = d ^ t2; t4 = a ^ t3; t5 = d | t1; t6 = t4 & t5; x1 = b ^ t6; t8 = t3 | x1; t9 = a & t4; x3 = t1 ^ t9; t11 = t4 ^ t8; t12 = x3 & t11; x2 = t3 ^ t12; x0 = (~t11) ^ (x3 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 8
		x0 ^= k[32];	x1 ^= k[33];	x2 ^= k[34];	x3 ^= k[35];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t3 = c ^ t1; t4 = b ^ t3; x3 = (a & d) ^ t4; t7 = a ^ (b & t1); x2 = t4 ^ (c | t7); t12 = x3 & (t3 ^ t7); x1 = (~t3) ^ t12; x0 = t12 ^ (~t7);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 9
		x0 ^= k[36];	x1 ^= k[37];	x2 ^= k[38];	x3 ^= k[39];
		a = x0; b = x1; c = x2; d = x3;
		t2 = b ^ (~a);		t5 = c ^ (a | t2);		x2 = d ^ t5;		t7 = b ^ (d | t2);		t8 = t2 ^ x2;		x3 = t8 ^ (t5 & t7); 		t11 = t5 ^ t7; 		x1 = x3 ^ t11; 		x0 = t5 ^ (t8 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);

		// Give difference
		ciphertext0 = x0; ciphertext1 = x1; ciphertext2 = x2; ciphertext3 = x3;
		x0 = plaintext0 ^ 0x20010000;
		x1 = plaintext1 ^ 0x00010000;
		x2 = plaintext2 ^ 0x00000000;
		x3 = plaintext3 ^ 0x20010000;
		//20010000000000000001000020010000

		// Round 4
		x0 ^= k[16];	x1 ^= k[17];	x2 ^= k[18];	x3 ^= k[19];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		round = 5;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 6
		round = 6;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 7
		round = 7;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = b ^ c; t2 = c & t1; t3 = d ^ t2; t4 = a ^ t3; t5 = d | t1; t6 = t4 & t5; x1 = b ^ t6; t8 = t3 | x1; t9 = a & t4; x3 = t1 ^ t9; t11 = t4 ^ t8; t12 = x3 & t11; x2 = t3 ^ t12; x0 = (~t11) ^ (x3 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 8
		x0 ^= k[32];	x1 ^= k[33];	x2 ^= k[34];	x3 ^= k[35];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t3 = c ^ t1; t4 = b ^ t3; x3 = (a & d) ^ t4; t7 = a ^ (b & t1); x2 = t4 ^ (c | t7); t12 = x3 & (t3 ^ t7); x1 = (~t3) ^ t12; x0 = t12 ^ (~t7);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 9
		x0 ^= k[36];	x1 ^= k[37];	x2 ^= k[38];	x3 ^= k[39];
		a = x0; b = x1; c = x2; d = x3;
		t2 = b ^ (~a);		t5 = c ^ (a | t2);		x2 = d ^ t5;		t7 = b ^ (d | t2);		t8 = t2 ^ x2;		x3 = t8 ^ (t5 & t7); 		t11 = t5 ^ t7; 		x1 = x3 ^ t11; 		x0 = t5 ^ (t8 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);


		x0 ^= ciphertext0;
		x1 ^= ciphertext1;
		x2 ^= ciphertext2;
		x3 ^= ciphertext3;

		t0 = x0 & 0x20108480;
		t1 = x1 & 0x0302842c;
		t2 = x2 & 0x01004000;
		t3 = x3 & 0x00008420;
		// 
		// corrected: 0010002600010020b4000b0010a02200 => 00008420010040000302842c20108480

		t0 ^= t0 >> 1;
		t0 ^= t0 >> 2;
		t0 = (t0 & 0x11111111U) * 0x11111111U;
		t0 = (t0 >> 28) & 1;

		t1 ^= t1 >> 1;
		t1 ^= t1 >> 2;
		t1 = (t1 & 0x11111111U) * 0x11111111U;
		t1 = (t1 >> 28) & 1;

		t2 ^= t2 >> 1;
		t2 ^= t2 >> 2;
		t2 = (t2 & 0x11111111U) * 0x11111111U;
		t2 = (t2 >> 28) & 1;

		t3 ^= t3 >> 1;
		t3 ^= t3 >> 2;
		t3 = (t3 & 0x11111111U) * 0x11111111U;
		t3 = (t3 >> 28) & 1;

		t0 = t0 ^ t1 ^ t2 ^ t3;

		if (t0 == 0) counter++;

		// Change plaintext
		plaintext0 = x0; plaintext1 = x1; plaintext2 = x2; plaintext3 = x3;
	}
	hit[blockIdx.x * blockDim.x + threadIdx.x] = counter;
}
__global__ void Serpent7round_PZWD24_middle4(bit32 p0, bit32 p1, bit32 p2, bit32 p3, bit64* hit, bit32* k_d) {
	__shared__ bit32  k[132];
	bit32 x0, x1, x2, x3, plaintext0 = p0, plaintext1 = p1, plaintext2 = p2, plaintext3 = p3 + blockIdx.x * blockDim.x + threadIdx.x;
	bit32 a, b, c, d;
	bit32 ciphertext0, ciphertext1, ciphertext2, ciphertext3;
	bit32 t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	bit32 round;
	bit64 counter = 0;
	if (threadIdx.x < 132) k[threadIdx.x] = k_d[threadIdx.x];
	__syncthreads();
	for (int i = 0; i < loop_d; i++) {
		x0 = plaintext0; x1 = plaintext1; x2 = plaintext2; x3 = plaintext3;
		// Round 2
		x0 ^= k[8];	x1 ^= k[9];	x2 ^= k[19];	x3 ^= k[11];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 3
		x0 ^= k[12];	x1 ^= k[13];	x2 ^= k[14];	x3 ^= k[15];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 4
		x0 ^= k[16];	x1 ^= k[17];	x2 ^= k[18];	x3 ^= k[19];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		x0 ^= k[20];	x1 ^= k[21];	x2 ^= k[22];	x3 ^= k[23];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);

		// Give difference
		ciphertext0 = x0; ciphertext1 = x1; ciphertext2 = x2; ciphertext3 = x3;
		x0 = plaintext0 ^ 0x00000000;
		x1 = plaintext1 ^ 0x00000000;
		x2 = plaintext2 ^ 0x00010000;
		x3 = plaintext3 ^ 0x00000000;
		//00000000000100000000000000000000

		// Round 2
		x0 ^= k[8];	x1 ^= k[9];	x2 ^= k[19];	x3 ^= k[11];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 3
		x0 ^= k[12];	x1 ^= k[13];	x2 ^= k[14];	x3 ^= k[15];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 4
		x0 ^= k[16];	x1 ^= k[17];	x2 ^= k[18];	x3 ^= k[19];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		x0 ^= k[20];	x1 ^= k[21];	x2 ^= k[22];	x3 ^= k[23];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);


		x0 ^= ciphertext0;
		x1 ^= ciphertext1;
		x2 ^= ciphertext2;
		x3 ^= ciphertext3;

		t0 = x0 & 0x00000000;
		t1 = x1 & 0x00000000;
		t2 = x2 & 0x00000000;
		t3 = x3 & 0x00000100;
		//00000100000000000000000000000000


		t0 ^= t0 >> 1;
		t0 ^= t0 >> 2;
		t0 = (t0 & 0x11111111U) * 0x11111111U;
		t0 = (t0 >> 28) & 1;

		t1 ^= t1 >> 1;
		t1 ^= t1 >> 2;
		t1 = (t1 & 0x11111111U) * 0x11111111U;
		t1 = (t1 >> 28) & 1;

		t2 ^= t2 >> 1;
		t2 ^= t2 >> 2;
		t2 = (t2 & 0x11111111U) * 0x11111111U;
		t2 = (t2 >> 28) & 1;

		t3 ^= t3 >> 1;
		t3 ^= t3 >> 2;
		t3 = (t3 & 0x11111111U) * 0x11111111U;
		t3 = (t3 >> 28) & 1;

		t0 = t0 ^ t1 ^ t2 ^ t3;

		if (t0 == 0) counter++;

		// Change plaintext
		plaintext0 = x0; plaintext1 = x1; plaintext2 = x2; plaintext3 = x3;
	}
	hit[blockIdx.x * blockDim.x + threadIdx.x] = counter;
}
__global__ void Serpent7round_PZWD24_middle5(bit32 p0, bit32 p1, bit32 p2, bit32 p3, bit64* hit, bit32* k_d) {
	__shared__ bit32  k[132];
	bit32 x0, x1, x2, x3, plaintext0 = p0, plaintext1 = p1, plaintext2 = p2, plaintext3 = p3 + blockIdx.x * blockDim.x + threadIdx.x;
	bit32 a, b, c, d;
	bit32 ciphertext0, ciphertext1, ciphertext2, ciphertext3;
	bit32 t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	bit32 round;
	bit64 counter = 0;
	if (threadIdx.x < 132) k[threadIdx.x] = k_d[threadIdx.x];
	__syncthreads();
	for (int i = 0; i < loop_d; i++) {
		x0 = plaintext0; x1 = plaintext1; x2 = plaintext2; x3 = plaintext3;
		// Round 2
		x0 ^= k[8];	x1 ^= k[9];	x2 ^= k[19];	x3 ^= k[11];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 3
		x0 ^= k[12];	x1 ^= k[13];	x2 ^= k[14];	x3 ^= k[15];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 4
		x0 ^= k[16];	x1 ^= k[17];	x2 ^= k[18];	x3 ^= k[19];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		x0 ^= k[20];	x1 ^= k[21];	x2 ^= k[22];	x3 ^= k[23];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 6
		x0 ^= k[24];	x1 ^= k[25];	x2 ^= k[26];	x3 ^= k[27];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);

		// Give difference
		ciphertext0 = x0; ciphertext1 = x1; ciphertext2 = x2; ciphertext3 = x3;
		x0 = plaintext0 ^ 0x00000000;
		x1 = plaintext1 ^ 0x00000000;
		x2 = plaintext2 ^ 0x00010000;
		x3 = plaintext3 ^ 0x00000000;
		//00000000000100000000000000000000

		// Round 2
		x0 ^= k[8];	x1 ^= k[9];	x2 ^= k[19];	x3 ^= k[11];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 3
		x0 ^= k[12];	x1 ^= k[13];	x2 ^= k[14];	x3 ^= k[15];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 4
		x0 ^= k[16];	x1 ^= k[17];	x2 ^= k[18];	x3 ^= k[19];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		x0 ^= k[20];	x1 ^= k[21];	x2 ^= k[22];	x3 ^= k[23];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 6
		x0 ^= k[24];	x1 ^= k[25];	x2 ^= k[26];	x3 ^= k[27];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);


		x0 ^= ciphertext0;
		x1 ^= ciphertext1;
		x2 ^= ciphertext2;
		x3 ^= ciphertext3;

		t0 = x0 & 0x04000000;
		t1 = x1 & 0x00200000;
		t2 = x2 & 0x00000000;
		t3 = x3 & 0x00200000;
		//00200000000000000020000004000000

		t0 ^= t0 >> 1;
		t0 ^= t0 >> 2;
		t0 = (t0 & 0x11111111U) * 0x11111111U;
		t0 = (t0 >> 28) & 1;

		t1 ^= t1 >> 1;
		t1 ^= t1 >> 2;
		t1 = (t1 & 0x11111111U) * 0x11111111U;
		t1 = (t1 >> 28) & 1;

		t2 ^= t2 >> 1;
		t2 ^= t2 >> 2;
		t2 = (t2 & 0x11111111U) * 0x11111111U;
		t2 = (t2 >> 28) & 1;

		t3 ^= t3 >> 1;
		t3 ^= t3 >> 2;
		t3 = (t3 & 0x11111111U) * 0x11111111U;
		t3 = (t3 >> 28) & 1;

		t0 = t0 ^ t1 ^ t2 ^ t3;

		if (t0 == 0) counter++;

		// Change plaintext
		plaintext0 = x0; plaintext1 = x1; plaintext2 = x2; plaintext3 = x3;
	}
	hit[blockIdx.x * blockDim.x + threadIdx.x] = counter;
}
__global__ void Serpent7round_PZWD24_middle6(bit32 p0, bit32 p1, bit32 p2, bit32 p3, bit64* hit, bit32* k_d) {
	__shared__ bit32  k[132];
	bit32 x0, x1, x2, x3, plaintext0 = p0, plaintext1 = p1, plaintext2 = p2, plaintext3 = p3 + blockIdx.x * blockDim.x + threadIdx.x;
	bit32 a, b, c, d;
	bit32 ciphertext0, ciphertext1, ciphertext2, ciphertext3;
	bit32 t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	bit32 round;
	bit64 counter = 0;
	if (threadIdx.x < 132) k[threadIdx.x] = k_d[threadIdx.x];
	__syncthreads();
	for (int i = 0; i < loop_d; i++) {
		x0 = plaintext0; x1 = plaintext1; x2 = plaintext2; x3 = plaintext3;
		// Round 2
		x0 ^= k[8];	x1 ^= k[9];	x2 ^= k[19];	x3 ^= k[11];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 3
		x0 ^= k[12];	x1 ^= k[13];	x2 ^= k[14];	x3 ^= k[15];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 4
		x0 ^= k[16];	x1 ^= k[17];	x2 ^= k[18];	x3 ^= k[19];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		x0 ^= k[20];	x1 ^= k[21];	x2 ^= k[22];	x3 ^= k[23];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 6
		x0 ^= k[24];	x1 ^= k[25];	x2 ^= k[26];	x3 ^= k[27];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 7
		x0 ^= k[28];	x1 ^= k[29];	x2 ^= k[30];	x3 ^= k[31];
		a = x0; b = x1; c = x2; d = x3;
		t1 = b ^ c; t2 = c & t1; t3 = d ^ t2; t4 = a ^ t3; t5 = d | t1; t6 = t4 & t5; x1 = b ^ t6; t8 = t3 | x1; t9 = a & t4; x3 = t1 ^ t9; t11 = t4 ^ t8; t12 = x3 & t11; x2 = t3 ^ t12; x0 = (~t11) ^ (x3 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);

		// Give difference
		ciphertext0 = x0; ciphertext1 = x1; ciphertext2 = x2; ciphertext3 = x3;
		x0 = plaintext0 ^ 0x00000000;
		x1 = plaintext1 ^ 0x00000000;
		x2 = plaintext2 ^ 0x00010000;
		x3 = plaintext3 ^ 0x00000000;
		//00000000000100000000000000000000

		// Round 2
		x0 ^= k[8];	x1 ^= k[9];	x2 ^= k[19];	x3 ^= k[11];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 3
		x0 ^= k[12];	x1 ^= k[13];	x2 ^= k[14];	x3 ^= k[15];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 4
		x0 ^= k[16];	x1 ^= k[17];	x2 ^= k[18];	x3 ^= k[19];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		x0 ^= k[20];	x1 ^= k[21];	x2 ^= k[22];	x3 ^= k[23];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 6
		x0 ^= k[24];	x1 ^= k[25];	x2 ^= k[26];	x3 ^= k[27];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 7
		x0 ^= k[28];	x1 ^= k[29];	x2 ^= k[30];	x3 ^= k[31];
		a = x0; b = x1; c = x2; d = x3;
		t1 = b ^ c; t2 = c & t1; t3 = d ^ t2; t4 = a ^ t3; t5 = d | t1; t6 = t4 & t5; x1 = b ^ t6; t8 = t3 | x1; t9 = a & t4; x3 = t1 ^ t9; t11 = t4 ^ t8; t12 = x3 & t11; x2 = t3 ^ t12; x0 = (~t11) ^ (x3 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);


		x0 ^= ciphertext0;
		x1 ^= ciphertext1;
		x2 ^= ciphertext2;
		x3 ^= ciphertext3;

		t0 = x0 & 0x84001080;
		t1 = x1 & 0x0c284084;
		t2 = x2 & 0x00090800;
		t3 = x3 & 0x20000084;
		//20000084000908000c28408484001080

		t0 ^= t0 >> 1;
		t0 ^= t0 >> 2;
		t0 = (t0 & 0x11111111U) * 0x11111111U;
		t0 = (t0 >> 28) & 1;

		t1 ^= t1 >> 1;
		t1 ^= t1 >> 2;
		t1 = (t1 & 0x11111111U) * 0x11111111U;
		t1 = (t1 >> 28) & 1;

		t2 ^= t2 >> 1;
		t2 ^= t2 >> 2;
		t2 = (t2 & 0x11111111U) * 0x11111111U;
		t2 = (t2 >> 28) & 1;

		t3 ^= t3 >> 1;
		t3 ^= t3 >> 2;
		t3 = (t3 & 0x11111111U) * 0x11111111U;
		t3 = (t3 >> 28) & 1;

		t0 = t0 ^ t1 ^ t2 ^ t3;

		if (t0 == 0) counter++;

		// Change plaintext
		plaintext0 = x0; plaintext1 = x1; plaintext2 = x2; plaintext3 = x3;
	}
	hit[blockIdx.x * blockDim.x + threadIdx.x] = counter;
}
__global__ void Serpent_benchmark(bit32 p0, bit32 p1, bit32 p2, bit32 p3, bit64* hit, bit32* k_d) {
	__shared__ bit32  k[132];
	bit32 x0, x1, x2, x3, plaintext0 = p0, plaintext1 = p1, plaintext2 = p2, plaintext3 = p3 + blockIdx.x * blockDim.x + threadIdx.x;
	bit32 ciphertext0;
	bit32 a, b, c, d;
	bit32 t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	bit64 counter = 0;
	if (threadIdx.x < 132) k[threadIdx.x] = k_d[threadIdx.x];
	__syncthreads();
	for (int i = 0; i < loop_d; i++) {
		x0 = plaintext0; x1 = plaintext1; x2 = plaintext2; x3 = plaintext3;
//#pragma unroll
/*		for (int i = 0; i < 32; i++) {
			x0 ^= k[i*4];	x1 ^= k[i * 4 + 1];	x2 ^= k[i * 4 + 2]; x3 ^= k[i * 4 + 3];
			a = x0; b = x1; c = x2; d = x3;
			t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
			x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		}*/

		// Round 0
		x0 ^= k[0];	x1 ^= k[1];	x2 ^= k[2];	x3 ^= k[3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t3 = c ^ t1; t4 = b ^ t3; x3 = (a & d) ^ t4; t7 = a ^ (b & t1); x2 = t4 ^ (c | t7); t12 = x3 & (t3 ^ t7); x1 = (~t3) ^ t12; x0 = t12 ^ (~t7);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 1
		x0 ^= k[4];	x1 ^= k[5];	x2 ^= k[6];	x3 ^= k[7];
		a = x0; b = x1; c = x2; d = x3;
		t2 = b ^ (~a);		t5 = c ^ (a | t2);		x2 = d ^ t5;		t7 = b ^ (d | t2);		t8 = t2 ^ x2;		x3 = t8 ^ (t5 & t7); 		t11 = t5 ^ t7; 		x1 = x3 ^ t11; 		x0 = t5 ^ (t8 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 2
		x0 ^= k[8];	x1 ^= k[9];	x2 ^= k[19];	x3 ^= k[11];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 3
		x0 ^= k[12];	x1 ^= k[13];	x2 ^= k[14];	x3 ^= k[15];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 4
		x0 ^= k[16];	x1 ^= k[17];	x2 ^= k[18];	x3 ^= k[19];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		x0 ^= k[20];	x1 ^= k[21];	x2 ^= k[22];	x3 ^= k[23];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 6
		x0 ^= k[24];	x1 ^= k[25];	x2 ^= k[26];	x3 ^= k[27];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 7
		x0 ^= k[28];	x1 ^= k[29];	x2 ^= k[30];	x3 ^= k[31];
		a = x0; b = x1; c = x2; d = x3;
		t1 = b ^ c; t2 = c & t1; t3 = d ^ t2; t4 = a ^ t3; t5 = d | t1; t6 = t4 & t5; x1 = b ^ t6; t8 = t3 | x1; t9 = a & t4; x3 = t1 ^ t9; t11 = t4 ^ t8; t12 = x3 & t11; x2 = t3 ^ t12; x0 = (~t11) ^ (x3 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 8
		x0 ^= k[32];	x1 ^= k[33];	x2 ^= k[34];	x3 ^= k[35];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t3 = c ^ t1; t4 = b ^ t3; x3 = (a & d) ^ t4; t7 = a ^ (b & t1); x2 = t4 ^ (c | t7); t12 = x3 & (t3 ^ t7); x1 = (~t3) ^ t12; x0 = t12 ^ (~t7);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 9
		x0 ^= k[36];	x1 ^= k[37];	x2 ^= k[38];	x3 ^= k[39];
		a = x0; b = x1; c = x2; d = x3;
		t2 = b ^ (~a);		t5 = c ^ (a | t2);		x2 = d ^ t5;		t7 = b ^ (d | t2);		t8 = t2 ^ x2;		x3 = t8 ^ (t5 & t7); 		t11 = t5 ^ t7; 		x1 = x3 ^ t11; 		x0 = t5 ^ (t8 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 10
		x0 ^= k[40];	x1 ^= k[41];	x2 ^= k[42];	x3 ^= k[43];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 11
		x0 ^= k[44];	x1 ^= k[45];	x2 ^= k[46];	x3 ^= k[47];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 12
		x0 ^= k[48];	x1 ^= k[49];	x2 ^= k[50];	x3 ^= k[51];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 13
		x0 ^= k[52];	x1 ^= k[53];	x2 ^= k[54];	x3 ^= k[55];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 14
		x0 ^= k[56];	x1 ^= k[57];	x2 ^= k[58];	x3 ^= k[59];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 15
		x0 ^= k[60];	x1 ^= k[61];	x2 ^= k[62];	x3 ^= k[63];
		a = x0; b = x1; c = x2; d = x3;
		t1 = b ^ c; t2 = c & t1; t3 = d ^ t2; t4 = a ^ t3; t5 = d | t1; t6 = t4 & t5; x1 = b ^ t6; t8 = t3 | x1; t9 = a & t4; x3 = t1 ^ t9; t11 = t4 ^ t8; t12 = x3 & t11; x2 = t3 ^ t12; x0 = (~t11) ^ (x3 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);

		// Round 16
		x0 ^= k[64];	x1 ^= k[65];	x2 ^= k[66];	x3 ^= k[67];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d;		t3 = c ^ t1;		t4 = b ^ t3;		x3 = (a & d) ^ t4;		t7 = a ^ (b & t1);		x2 = t4 ^ (c | t7);		t12 = x3 & (t3 ^ t7);		x1 = (~t3) ^ t12;		x0 = t12 ^ (~t7);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 17
		x0 ^= k[68];	x1 ^= k[69];	x2 ^= k[70];	x3 ^= k[71];
		a = x0; b = x1; c = x2; d = x3;
		t2 = b ^ (~a);		t5 = c ^ (a | t2);		x2 = d ^ t5;		t7 = b ^ (d | t2);		t8 = t2 ^ x2;		x3 = t8 ^ (t5 & t7); 		t11 = t5 ^ t7; 		x1 = x3 ^ t11; 		x0 = t5 ^ (t8 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 18
		x0 ^= k[72];	x1 ^= k[73];	x2 ^= k[74];	x3 ^= k[75];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 19
		x0 ^= k[76];	x1 ^= k[77];	x2 ^= k[78];	x3 ^= k[79];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 20
		x0 ^= k[80];	x1 ^= k[81];	x2 ^= k[82];	x3 ^= k[83];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 21
		x0 ^= k[84];	x1 ^= k[85];	x2 ^= k[86];	x3 ^= k[87];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 22
		x0 ^= k[88];	x1 ^= k[89];	x2 ^= k[90];	x3 ^= k[91];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 23
		x0 ^= k[92];	x1 ^= k[93];	x2 ^= k[94];	x3 ^= k[95];
		a = x0; b = x1; c = x2; d = x3;
		t1 = b ^ c; t2 = c & t1; t3 = d ^ t2; t4 = a ^ t3; t5 = d | t1; t6 = t4 & t5; x1 = b ^ t6; t8 = t3 | x1; t9 = a & t4; x3 = t1 ^ t9; t11 = t4 ^ t8; t12 = x3 & t11; x2 = t3 ^ t12; x0 = (~t11) ^ (x3 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 24
		x0 ^= k[96];	x1 ^= k[97];	x2 ^= k[98];	x3 ^= k[99];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t3 = c ^ t1; t4 = b ^ t3; x3 = (a & d) ^ t4; t7 = a ^ (b & t1); x2 = t4 ^ (c | t7); t12 = x3 & (t3 ^ t7); x1 = (~t3) ^ t12; x0 = t12 ^ (~t7);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 25
		x0 ^= k[100];	x1 ^= k[101];	x2 ^= k[102];	x3 ^= k[103];
		a = x0; b = x1; c = x2; d = x3;
		t2 = b ^ (~a);		t5 = c ^ (a | t2);		x2 = d ^ t5;		t7 = b ^ (d | t2);		t8 = t2 ^ x2;		x3 = t8 ^ (t5 & t7); 		t11 = t5 ^ t7; 		x1 = x3 ^ t11; 		x0 = t5 ^ (t8 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 26
		x0 ^= k[104];	x1 ^= k[105];	x2 ^= k[106];	x3 ^= k[107];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 27
		x0 ^= k[108];	x1 ^= k[109];	x2 ^= k[110];	x3 ^= k[111];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 28
		x0 ^= k[112];	x1 ^= k[113];	x2 ^= k[114];	x3 ^= k[115];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 29
		x0 ^= k[116];	x1 ^= k[117];	x2 ^= k[118];	x3 ^= k[119];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 30
		x0 ^= k[120];	x1 ^= k[121];	x2 ^= k[122];	x3 ^= k[123];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 31
		x0 ^= k[124];	x1 ^= k[125];	x2 ^= k[126];	x3 ^= k[127];
		a = x0; b = x1; c = x2; d = x3;
		t1 = b ^ c; t2 = c & t1; t3 = d ^ t2; t4 = a ^ t3; t5 = d | t1; t6 = t4 & t5; x1 = b ^ t6; t8 = t3 | x1; t9 = a & t4; x3 = t1 ^ t9; t11 = t4 ^ t8; t12 = x3 & t11; x2 = t3 ^ t12; x0 = (~t11) ^ (x3 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 31
		x0 ^= k[128];	x1 ^= k[129];	x2 ^= k[130];	x3 ^= k[131];

		if (x0 == 0)
			if (x1 == 0)
				if (x2 == 0)
					if (x3 == 0)
						printf("Hello world\n");

		plaintext0 = x0; plaintext1 = x1; plaintext2 = x2; plaintext3 = x3;
	}
	hit[blockIdx.x * blockDim.x + threadIdx.x] = counter;
}
__global__ void Serpent5round_eprint(bit32 p0, bit32 p1, bit32 p2, bit32 p3, bit64* hit, bit32* k_d) {
	__shared__ bit32  k[132];
	bit32 x0, x1, x2, x3, plaintext0 = p0, plaintext1 = p1, plaintext2 = p2, plaintext3 = p3 + blockIdx.x * blockDim.x + threadIdx.x;
	bit32 ciphertext0, ciphertext1, ciphertext2, ciphertext3;
	bit32 a, b, c, d;
	bit32 t0, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10, t11, t12;
	bit32 round;
	bit64 counter = 0;
	if (threadIdx.x < 132) k[threadIdx.x] = k_d[threadIdx.x];
	__syncthreads();
	for (int i = 0; i < loop_d; i++) {
		x0 = plaintext0; x1 = plaintext1; x2 = plaintext2; x3 = plaintext3;
		// Round 2
		round = 2;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 3
		round = 3;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 4
		round = 4;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		round = 5;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 6
		round = 6;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		// no linear transformation for the last round

			// Last round linear transformation
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);

		ciphertext0 = x0;
		ciphertext1 = x1;
		ciphertext2 = x2;
		ciphertext3 = x3;
		// Give difference	(eprint)	

		x0 = plaintext0 ^ (0x01000000);
		x1 = plaintext1 ^ (0x01000000);
		x2 = plaintext2 ^ (0x09000000);
		x3 = plaintext3 ^ (0x00000000);

		// Round 2
		round = 2;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = b ^ d; t3 = c & t1; x0 = t2 ^ t3; t5 = c ^ t1; t6 = c ^ x0; t7 = b & t6; x3 = t5 ^ t7; x2 = a ^ ((d | t7) & (x0 | t5)); x1 = (t2 ^ x3) ^ (x2 ^ (d | t1));
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 3
		round = 3;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ b; t2 = a & c; t3 = a | d; t4 = c ^ d; t5 = t1 & t3; t6 = t2 | t5; x2 = t4 ^ t6; t8 = b ^ t3; t9 = t6 ^ t8; t10 = t4 & t9; x0 = t1 ^ t10; t12 = x2 & x0; x1 = t9 ^ t12; x3 = (b | d) ^ (t4 ^ t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 4
		round = 4;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = a ^ d; t2 = d & t1; t3 = c ^ t2; t4 = b | t3; x3 = t1 ^ t4; t6 = ~b; t7 = t1 | t6; x0 = t3 ^ t7; t9 = a & x0; t10 = t1 ^ t6; t11 = t4 & t10; x2 = t9 ^ t11; x1 = (a ^ t3) ^ (t10 & x2);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 5
		round = 5;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ b; t3 = a ^ d; t4 = c ^ t1; t5 = t2 | t3; x0 = t4 ^ t5; t7 = d & x0; t8 = t2 ^ x0; x1 = t7 ^ t8; t10 = t1 | x0; t11 = t2 | t7; t12 = t3 ^ t10; x2 = t11 ^ t12; x3 = (b ^ t7) ^ (x1 & t12);
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);
		// Round 6
		round = 6;
		x0 = x0 ^ k[4 * round];	x1 = x1 ^ k[4 * round + 1];	x2 = x2 ^ k[4 * round + 2];	x3 = x3 ^ k[4 * round + 3];
		a = x0; b = x1; c = x2; d = x3;
		t1 = ~a; t2 = a ^ d; t3 = b ^ t2; t4 = t1 | t2; t5 = c ^ t4; x1 = b ^ t5; t7 = t2 | x1; t8 = d ^ t7; t9 = t5 & t8; x2 = t3 ^ t9; t11 = t5 ^ t8; x0 = x2 ^ t11; x3 = (~t5) ^ (t3 & t11);
		// no linear transformation for the last round

			// Last round linear transformation
		x0 = (x0 << 13) | x0 >> 19;	x2 = (x2 << 3) | (x2 >> 29); x1 = x1 ^ x0 ^ x2;	x3 = x3 ^ x2 ^ (x0 << 3); x1 = (x1 << 1) | (x1 >> 31); x3 = (x3 << 7) | (x3 >> 25);	x0 = x0 ^ x1 ^ x3; x2 = x2 ^ x3 ^ (x1 << 7); x0 = (x0 << 5) | (x0 >> 27); x2 = (x2 << 22) | (x2 >> 10);

		// Check Difference
		x0 ^= ciphertext0;
		x1 ^= ciphertext1;
		x2 ^= ciphertext2;
		x3 ^= ciphertext3;

		// Masked bits (eprint)
/*		t0 = x0 & 0x84001080;
		t1 = x1 & 0x0c284084;
		t2 = x2 & 0x00090800;
		t3 = x3 & 0x20000084;*/

		// without the final transformation
/*		t0 = x0 & 0x00000004;
		t1 = x1 & 0x00000004;
		t2 = x2 & 0x00000000;
		t3 = x3 & 0x00000004;*/

		t0 = x0 & 0x42004000;
		t1 = x1 & 0x02100228;
		t2 = x2 & 0x00000005;
		t3 = x3 & 0x02181600;



		// Our best mask
/*		t0 = x0 & 0x02100000;
		t1 = x1 & 0x00000000;
		t2 = x2 & 0x00000000;
		t3 = x3 & 0x00000000;*/

		t0 ^= t0 >> 1;
		t0 ^= t0 >> 2;
		t0 = (t0 & 0x11111111U) * 0x11111111U;
		t0 = (t0 >> 28) & 1;

		t1 ^= t1 >> 1;
		t1 ^= t1 >> 2;
		t1 = (t1 & 0x11111111U) * 0x11111111U;
		t1 = (t1 >> 28) & 1;

		t2 ^= t2 >> 1;
		t2 ^= t2 >> 2;
		t2 = (t2 & 0x11111111U) * 0x11111111U;
		t2 = (t2 >> 28) & 1;

		t3 ^= t3 >> 1;
		t3 ^= t3 >> 2;
		t3 = (t3 & 0x11111111U) * 0x11111111U;
		t3 = (t3 >> 28) & 1;

		t0 = t0 ^ t1 ^ t2 ^ t3;

		if (t0 == 0) counter++;

		// Change plaintext
		plaintext0 = x0; plaintext1 = x1; plaintext2 = x2; plaintext3 = x3;
	}
	hit[blockIdx.x * blockDim.x + threadIdx.x] = counter;
}

int main() {	
	FILE* fp3, * fp, * fp2, * fp4, * fp5;
	bit64* hit_d;
	bit64* hit_test;
	__int64 temp, hit, bias = 0, cumulative_bias[100] = { 0 }, average_bias;
	bit32* k_d, keylist[100][4];
	hit_test = (bit64*)calloc(totalthreads, sizeof(bit64));
	int t, j, i;
	cudaSetDeviceFlags(cudaDeviceScheduleBlockingSync);
	get_user_inputs();
	cudaSetDevice(devicenumber);
	cudaMalloc((void**)&hit_d, totalthreads * sizeof(bit64));
	cudaMalloc((void**)&k_d, 132 * sizeof(bit32));
	if (numberofrounds == 3) fp = fopen("3round_cumulative.txt", "ab");
	else if (numberofrounds == 4) fp = fopen("4round_cumulative.txt", "ab");
	else if (numberofrounds == 5) fp = fopen("5round_cumulative.txt", "ab");
	else if (numberofrounds == 66) fp = fopen("66round_cumulative.txt", "ab");
	else if (numberofrounds == 74) fp = fopen("74round_cumulative.txt", "ab");
	else if (numberofrounds == 75) fp = fopen("75round_cumulative.txt", "ab");
	else if (numberofrounds == 76) fp = fopen("76round_cumulative.txt", "ab");
	else if (numberofrounds == 32) fp = fopen("32round_benchmark.txt", "ab");
	else exit(1);
	fprintf(fp, "Pairs: %I64d\n", blocksize * threadsize * loop * pairs);
	plaintext[0] = 0xc7675a6e; plaintext[1] = 0xe3d21628; plaintext[2] = 0x9f090c02; plaintext[3] = 0x5eeaf89d;
	fp3 = fopen("keys.txt", "r");
	for (i = 0; i < 100; i++) { fscanf(fp3, "%x", &keylist[i][0]);	fscanf(fp3, "%x", &keylist[i][1]);	fscanf(fp3, "%x", &keylist[i][2]);	fscanf(fp3, "%x", &keylist[i][3]); }
	fclose(fp3);
	if (startingpoint) {
		fp4 = fopen("startingpoint.txt", "r");	for (i = 0; i < experiment; i++) { fscanf(fp4, "%I64d", &cumulative_bias[i]); printf("%d: %I64d\n", i, cumulative_bias[i]); }	fclose(fp4);
		for (i = 0; i < startingpoint; i++) { plaintext[0]++; plaintext[1]++; plaintext[2]++; plaintext[3]++; }
	}
	for (j = 0; j < pairs; j++) {
		StartCounter();
		for (t = 0; t < experiment; t++) {
			hit = 0;
			key[0] = keylist[t][0]; key[1] = keylist[t][1]; key[2] = keylist[t][2]; key[3] = keylist[t][3];
			key_schedule();
			cudaMemcpy(k_d, k, 132 * sizeof(bit32), cudaMemcpyHostToDevice);
			cudaMemcpy(hit_d, hit_test, totalthreads * sizeof(bit64), cudaMemcpyHostToDevice);
			
			if (numberofrounds == 3) Serpent3round_PZWD24 << <blocksize, threadsize >> > (plaintext[0], plaintext[1], plaintext[2], plaintext[3], hit_d, k_d);
			else if (numberofrounds == 4) Serpent4round << <blocksize, threadsize >> > (plaintext[0], plaintext[1], plaintext[2], plaintext[3], hit_d, k_d);
			else if (numberofrounds == 5) Serpent5round_eprint << <blocksize, threadsize >> > (plaintext[0], plaintext[1], plaintext[2], plaintext[3], hit_d, k_d);
			else if (numberofrounds == 66) Serpent6round_PZWD24_final << <blocksize, threadsize >> > (plaintext[0], plaintext[1], plaintext[2], plaintext[3], hit_d, k_d);
			else if (numberofrounds == 74) Serpent7round_PZWD24_middle4 << <blocksize, threadsize >> > (plaintext[0], plaintext[1], plaintext[2], plaintext[3], hit_d, k_d);
			else if (numberofrounds == 75) Serpent7round_PZWD24_middle5 << <blocksize, threadsize >> > (plaintext[0], plaintext[1], plaintext[2], plaintext[3], hit_d, k_d);
			else if (numberofrounds == 76) Serpent7round_PZWD24_middle6 << <blocksize, threadsize >> > (plaintext[0], plaintext[1], plaintext[2], plaintext[3], hit_d, k_d);
			else if (numberofrounds == 32) {
				float milliseconds = 0;
				cudaDeviceSynchronize();
				cudaEvent_t start, stop;	cudaEventCreate(&start);	cudaEventCreate(&stop);	cudaEventRecord(start);
				Serpent_benchmark << <blocksize, threadsize >> > (plaintext[0], plaintext[1], plaintext[2], plaintext[3], hit_d, k_d);
				cudaEventRecord(stop);	cudaEventSynchronize(stop);
				cudaEventElapsedTime(&milliseconds, start, stop);	printf("Time elapsed: %f milliseconds\n", milliseconds);
				exit(1);
			}
			cudaMemcpy(hit_test, hit_d, totalthreads * sizeof(bit64), cudaMemcpyDeviceToHost);
			for (i = 0; i < totalthreads; i++) hit += hit_test[i];
			cumulative_bias[t] += (hit - totalthreads*loop_d/2);
			//			printf("%2d: Key bias %I64d\n", t, cumulative_bias[t]);
			fprintf(fp, "%I64d\n", cumulative_bias[t]);
		}
		average_bias = 0;
		for (i = 0; i < experiment; i++) average_bias += cumulative_bias[i];
		average_bias /= experiment;
		printf("%5d: Time: %lf Average Bias %I64d\n", j + startingpoint, GetCounter(), average_bias);
		printf("\nTotal counter: %I64d\nDifference from the Expected Value: %I64d\nBias: 2^-%lf (For an experiment with 2^%lf data)\n", hit, average_bias, ((log(totalthreads * loop_d * (j+1))) / log(2)) - (log(abs(average_bias)) / log(2)), (log(totalthreads * loop_d * (j+1))) / log(2));
		fprintf(fp, "%5d: Time: %lf Average Bias %I64d\n", j + startingpoint, GetCounter(), average_bias);
		if (numberofrounds == 3) fp2 = fopen("3round_cumulative_bias.txt", "ab");
		else if (numberofrounds == 4) fp2 = fopen("4round_cumulative_bias.txt", "ab");
		else if (numberofrounds == 5) fp2 = fopen("5round_cumulative_bias.txt", "ab");
		else if (numberofrounds == 66) fp2 = fopen("66round_cumulative_bias.txt", "ab");
		else if (numberofrounds == 74) fp2 = fopen("74round_cumulative_bias.txt", "ab");
		else if (numberofrounds == 75) fp2 = fopen("75round_cumulative_bias.txt", "ab");
		else if (numberofrounds == 76) fp2 = fopen("76round_cumulative_bias.txt", "ab");
		fprintf(fp2, "%5d: Time: %lf Average Bias %I64d\n", j + startingpoint, GetCounter(), average_bias);
		fclose(fp2);
		fp5 = fopen("6round_cumulative_bias_shortened.txt", "ab");
		fprintf(fp5, "%I64d\n", average_bias);
		fclose(fp5);
		plaintext[0]++; plaintext[1]++; plaintext[2]++; plaintext[3]++;
	}
	cudaFree(hit_d); cudaFree(k_d);
	fclose(fp);
	free(hit_test);
	printf("%s\n", cudaGetErrorString(cudaGetLastError()));
	return 0;
}

