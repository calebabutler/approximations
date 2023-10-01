
// Includes
#include <stdio.h>
#include <stdint.h>
#include <math.h>

// Function prototypes
static double Floord(double);

static double SinStage1(double);
static double SinStage2(double);
static double SinStage3(double);

static double AtanStage1(double);
static double AtanStage2(double);
static double Atan(double);

static double Exp2Stage1(double);
static double Exp2(double);

static double Log2Stage1(double);
static double Log2(double);

// Global constants
static const double PI = 3.1415926535897932384626433832795028841971693993751058;
static const double DIV_2_PI = 1.0 / (2.0 * PI);
static const double INF = 1.0/0.0;
static const double DOUBLE_NAN = 0.0/0.0;

static const double SQRT2 = 1.4142135623730950488016887242096980785696718753769;

static const double LOG2E = 1.4426950408889634073599246810018921374266459541529;
static const double LOGE2 = 0.6931471805599453094172321214581765680755001343602;

double Floord(double x) {
	if (x >= 0)
		return (double) ((int) x);
	return (double) (((int) x) - 1);
}

/************
 * Math_Sin *
 ************/

double SinStage1(double x) {
	const double A[] = {
		.52359877559829885532,
		-.2392459620393377657e-1,
		.32795319441392666e-3,
		-.214071970654441e-5,
		.815113605169e-8,
		-.2020852964e-10,
	};

	double P = A[5];
	double x_2 = x * x;

	for (int i = 4; i >= 0; i--) {
		P *= x_2;
		P += A[i];
	}
	P *= x;
	return P;
}

double SinStage2(double x) {
	double sin_6 = SinStage1(x * 4.0);
	return sin_6 * (3.0 - 4.0 * sin_6 * sin_6);
}

double SinStage3(double x) {
	if (x < 0.25)
		return SinStage2(x);
	if (x < 0.5)
		return SinStage2(0.5 - x);
	if (x < 0.75)
		return -SinStage2(x - 0.5);
	return -SinStage2(1.0 - x);
}


double Math_Sin(double x) {
	double x_div_pi = x * DIV_2_PI;
	return SinStage3(x_div_pi - Floord(x_div_pi));
}

/************
 * Math_Cos *
 ************/

double Math_Cos(double x) {
	double x_div_pi_shifted = x * DIV_2_PI + 0.25;
	return SinStage3(x_div_pi_shifted - Floord(x_div_pi_shifted));
}

/**************
 * Math_Atan2 *
 **************/

double AtanStage1(double x) {
	const double A[] = {
		.99999999999969557,
		-.3333333333318,
		.1999999997276,
		-.14285702288,
		.11108719478,
		-.8870580341e-1,
	};

	double P = A[5];
	double x_2 = x * x;

	for (int i = 4; i >= 0; i--) {
		P *= x_2;
		P += A[i];
	}
	P *= x;
	return P;
}

double AtanStage2(double x) {
	const double X_i[] = {
		0.0,
		0.0984914033571642477671304050090839155018329620361328125,
		0.3033466836073424044428747947677038609981536865234375,
		0.53451113595079158269385288804187439382076263427734375,
		0.82067879082866024287312711749109439551830291748046875,
		1.218503525587976366040265929768793284893035888671875,
		1.8708684117893887854933154812897555530071258544921875,
		3.29655820893832096629694206058047711849212646484375,
		INF,
	};

	const double div_x_i[] = {
		0,
		0,
		5.02733949212584807497705696732737123966217041015625,
		2.41421356237309492343001693370752036571502685546875,
		1.496605762665489169904731170390732586383819580078125,
		1.0000000000000002220446049250313080847263336181640625,
		0.66817863791929898997778991542872972786426544189453125,
		0.414213562373095089963470627481001429259777069091796875,
		0.1989123673796580893391450217677629552781581878662109375,
	};

	const double div_x_i_2_plus_1[] = {
		0,
		0,
		26.2741423690881816810360760428011417388916015625,
		6.8284271247461898468600338674150407314300537109375,
		3.23982880884355051165357508580200374126434326171875,
		2.000000000000000444089209850062616169452667236328125,
		1.446462692171689656817079594475217163562774658203125,
		1.1715728752538099310953612075536511838436126708984375,
		1.0395661298965801488947136022034101188182830810546875,
	};

	int L = 0;
	int R = 8;

	while (R - L > 1) {
		int m = (L + R) / 2;
		if (X_i[m] <= x)
			L = m;
		else if (X_i[m] > x)
			R = m;
	}

	if (R <= 1)
		return AtanStage1(x);

	double t = div_x_i[R] - div_x_i_2_plus_1[R] / (div_x_i[R] + x);
	if (t >= 0)
		return (2 * R - 2) * PI / 32.0 + AtanStage1(t);

	return (2 * R - 2) * PI / 32.0 - AtanStage1(-t);
}

double Atan(double x) {
	if (x >= 0)
		return AtanStage2(x);
	return -AtanStage2(-x);
}

double Math_Atan2(double y, double x) {
	if (x > 0)
		return Atan(y / x);
	if (x < 0) {
		if (y >= 0)
			return Atan(y / x) + PI;
		return Atan(y / x) - PI;
	}
	if (y > 0)
		return PI / 2.0;
	if (y < 0)
		return -PI / 2.0;
	return DOUBLE_NAN;
}

/*************
 * Math_Exp2 *
 *************/

// EXPB 1067
double Exp2Stage1(double x) {
	const double A_P[] = {
		.1513906799054338915894328e4,
		.20202065651286927227886e2,
		.23093347753750233624e-1,
	};

	const double A_Q[] = {
		.4368211662727558498496814e4,
		.233184211427481623790295e3,
		1.0,
	};

	double x_2 = x * x;

	double P = A_P[2];
	for (int i = 1; i >= 0; i--) {
		P *= x_2;
		P += A_P[i];
	}
	P *= x;

	double Q = A_Q[2];
	for (int i = 1; i >= 0; i--) {
		Q *= x_2;
		Q += A_Q[i];
	}

	return (Q + P) / (Q - P);
}

double Exp2(double x) {
	int x_int = (int) x;

	if (x < 0)
		x_int--;

	if (x_int < -1022)
		return 0.0;
	if (x_int > 1023)
		return INF;

	union { double d; uint64_t i; } doi;
	doi.i = x_int + 1023;
	doi.i <<= 52;

	return doi.d * SQRT2 * Exp2Stage1(x - (double) x_int - 0.5);
}

/*************
 * Math_Log2 *
 *************/

// LOG2 2524
double Log2Stage1(double x) {
	const double A_P[] = {
		-.205466671951e1,
		-.88626599391e1,
		.610585199015e1,
		.481147460989e1,
	};

	const double A_Q[] = {
		.353553425277,
		.454517087629e1,
		.642784209029e1,
		1.0,
	};

	double P = A_P[3];
	for (int i = 2; i >= 0; i--) {
		P *= x;
		P += A_P[i];
	}

	double Q = A_Q[3];
	for (int i = 2; i >= 0; i--) {
		Q *= x;
		Q += A_Q[i];
	}

	return P / Q;
}

double Log2(double x) {
	if (x <= 0.0)
		return DOUBLE_NAN;

	union { double d; uint64_t i; } doi;
	doi.d = x;
	int integer_part = (doi.i >> 52);
	integer_part -= 1023;

	doi.i |= (((uint64_t) 1023) << 52);
	doi.i &= ~(((uint64_t) 1024) << 52);

	return integer_part + Log2Stage1(doi.d);
}

/************
 * Math_Exp *
 ************/

double Math_Exp(double x) {
	return Exp2(x * LOG2E);
}

/************
 * Math_Log *
 ************/

double Math_Log(double x) {
	return Log2(x) * LOGE2;
}

/*****************
 * main function *
 *****************/

int main(void) {
	double step = 20.0/1000000.0;

	for (double x = -10.0; x <= 10.0; x += step) {
#if defined(SIN)
		printf("%.60e %.60e\n", x, Math_Sin(x));
#elif defined(LIBM_SIN)
		printf("%.60e %.60e\n", x, sin(x));
#elif defined(COS)
		printf("%.60e %.60e\n", x, Math_Cos(x));
#elif defined(LIBM_COS)
		printf("%.60e %.60e\n", x, cos(x));
#elif defined(ATAN)
		printf("%.60e %.60e\n", x, Atan(x));
#elif defined(LIBM_ATAN)
		printf("%.60e %.60e\n", x, atan(x));
#elif defined(EXP)
		printf("%.60e %.60e\n", x, Math_Exp(x));
#elif defined(LIBM_EXP)
		printf("%.60e %.60e\n", x, exp(x));
#elif defined(LOG)
		printf("%.60e %.60e\n", x, Math_Log(x));
#elif defined(LIBM_LOG)
		printf("%.60e %.60e\n", x, log(x));
#endif
	}

	return 0;
}

