
/* Copyright (c) 2023 Caleb Butler
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

/* This code implements the math functions sine, cosine, arctangent, the
 * exponential function, and the logarithmic function. The code uses techniques
 * exclusively described in the book "Computer Approximations" by John Fraser
 * Hart (1st Edition). Each function approximates their associated math function
 * the same way:
 *
 *   1. First, the function uses properties of the associated math function to
 *      reduce the input range to a small finite interval,
 *
 *   2. Second, the function calculates a polynomial, rational, or similar
 *      function that approximates the associated math function on that small
 *      finite interval to the desired accuracy. These polynomial, rational, or
 *      similar functions were calculated by the authors of "Computer
 *      Approximations" using the Remez algorithm and exist in the book's
 *      appendix.
 */

/* Includes */
#include <stdio.h>

#if defined(LIBM_SIN) || defined(LIBM_COS) || defined(LIBM_ATAN) || defined(LIBM_EXP) || defined(LIBM_LOG)
#include <math.h>
#endif

/* Type definitions */
#if __MSC_VER
	typedef unsigned __int64 cc_uint64;
#elif __GNUC__
	typedef __UINT64_TYPE__ cc_uint64;
#else
	typedef unsigned long long cc_uint64;
#endif

/* Function prototypes */
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

static double Math_Sin(double);
static double Math_Cos(double);
static double Math_Atan2(double, double);
static double Math_Exp(double);
static double Math_Log(double);

/* Global constants */
static const double PI = 3.1415926535897932384626433832795028841971693993751058;
static const double DIV_2_PI = 1.0 / (2.0 * PI);
static const double INF = 1.0/0.0;
static const double DOUBLE_NAN = 0.0/0.0;

static const double SQRT2 = 1.4142135623730950488016887242096980785696718753769;

static const double LOG2E = 1.4426950408889634073599246810018921374266459541529;
static const double LOGE2 = 0.6931471805599453094172321214581765680755001343602;

/* Calculates the floor of a double.
 */
double Floord(double x) {
	if (x >= 0)
		return (double) ((int) x);
	return (double) (((int) x) - 1);
}

/************
 * Math_Sin *
 ************/

/* Calculates the 5th degree polynomial function SIN 2922 listed in the book's
 * appendix.
 *
 * Associated math function: sin(pi/6 * x)
 * Allowed input range: [0, 1]
 * Precision: 16.47
 */
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
	int i;

	for (i = 4; i >= 0; i--) {
		P *= x_2;
		P += A[i];
	}
	P *= x;
	return P;
}

/* Uses the property
 *   sin(x) = sin(x/3) * (3 - 4 * (sin(x/3))^2)
 * to reduce the input range of sin(x) to [0, pi/6].
 *
 * Associated math function: sin(2 * pi * x)
 * Allowed input range: [0, 0.25]
 */
double SinStage2(double x) {
	double sin_6 = SinStage1(x * 4.0);
	return sin_6 * (3.0 - 4.0 * sin_6 * sin_6);
}

/* Uses the properties of sine to reduce the input range from [0, 2*pi] to [0,
 * pi/2].
 *
 * Associated math function: sin(2 * pi * x)
 * Allowed input range: [0, 1]
 */
double SinStage3(double x) {
	if (x < 0.25)
		return SinStage2(x);
	if (x < 0.5)
		return SinStage2(0.5 - x);
	if (x < 0.75)
		return -SinStage2(x - 0.5);
	return -SinStage2(1.0 - x);
}

/* Since sine has a period of 2*pi, this function maps any real number to a
 * number from [0, 2*pi].
 *
 * Associated math function: sin(x)
 * Allowed input range: anything
 */
double Math_Sin(double x) {
	double x_div_pi = x * DIV_2_PI;
	return SinStage3(x_div_pi - Floord(x_div_pi));
}

/************
 * Math_Cos *
 ************/

/* This function works just like the above sine function, except it shifts the
 * input by pi/2, using the property cos(x) = sin(x + pi/2).
 *
 * Associated math function: cos(x)
 * Allowed input range: anything
 */
double Math_Cos(double x) {
	double x_div_pi_shifted = x * DIV_2_PI + 0.25;
	return SinStage3(x_div_pi_shifted - Floord(x_div_pi_shifted));
}

/**************
 * Math_Atan2 *
 **************/

/* Calculates the 5th degree polynomial ARCTN 4903 listed in the book's
 * appendix.
 *
 * Associated math function: arctan(x)
 * Allowed input range: [0, tan(pi/32)]
 * Precision: 16.52
 */
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
	int i;

	for (i = 4; i >= 0; i--) {
		P *= x_2;
		P += A[i];
	}
	P *= x;
	return P;
}

/* This function finds out in which partition the non-negative real number x
 * resides out of 8 partitions, which are precomputed. It then uses the
 * following law:
 *
 *   t = x_i^{-1} - (x_i^{-2} + 1)/(x_i^{-1} + x)
 *   arctan(x) = arctan(x_i) + arctan(t)
 *
 * where x_i = tan((2i - 2)*pi/32) and i is the partition number. The value of t
 * is guaranteed to be between [-tan(pi/32), tan(pi/32)].
 *
 * Associated math function: arctan(x)
 * Allowed input range: [0, infinity]
 */
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
	double t;

	while (R - L > 1) {
		int m = (L + R) / 2;
		if (X_i[m] <= x)
			L = m;
		else if (X_i[m] > x)
			R = m;
	}

	if (R <= 1)
		return AtanStage1(x);

	t = div_x_i[R] - div_x_i_2_plus_1[R] / (div_x_i[R] + x);
	if (t >= 0)
		return (2 * R - 2) * PI / 32.0 + AtanStage1(t);

	return (2 * R - 2) * PI / 32.0 - AtanStage1(-t);
}

/* Uses the property arctan(x) = -arctan(-x).
 *
 * Associated math function: arctan(x)
 * Allowed input range: anything
 */
double Atan(double x) {
	if (x >= 0)
		return AtanStage2(x);
	return -AtanStage2(-x);
}

/* Implements the function atan2 using Atan.
 *
 * Associated math function: atan2(y, x)
 * Allowed input range: anything
 */
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

/************
 * Math_Exp *
 ************/

/* Calculates the function EXPB 1067 listed in the book's appendix. It is of the
 * form
 *   (Q(x^2) + x*P(x^2)) / (Q(x^2) - x*P(x^2))
 *
 * Associated math function: 2^x
 * Allowed input range: [-1/2, 1/2]
 * Precision: 18.08
 */
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
	double P, Q;
	int i;

	P = A_P[2];
	for (i = 1; i >= 0; i--) {
		P *= x_2;
		P += A_P[i];
	}
	P *= x;

	Q = A_Q[2];
	for (i = 1; i >= 0; i--) {
		Q *= x_2;
		Q += A_Q[i];
	}

	return (Q + P) / (Q - P);
}

/* Reduces the range of 2^x to [-1/2, 1/2] by using the property
 *   2^x = 2^(integer value) * 2^(fractional part).
 * 2^(integer value) can be calculated by directly manipulating the bits of the
 * double-precision floating point representation.
 *
 * Associated math function: 2^x
 * Allowed input range: anything
 */
double Exp2(double x) {
	int x_int = (int) x;
	union { double d; cc_uint64 i; } doi;

	if (x < 0)
		x_int--;

	if (x_int < -1022)
		return 0.0;
	if (x_int > 1023)
		return INF;

	doi.i = x_int + 1023;
	doi.i <<= 52;

	return doi.d * SQRT2 * Exp2Stage1(x - (double) x_int - 0.5);
}

/* Uses the fact that
 *   exp(x) = 2^(x * log_2(e)).
 *
 * Associated math function: exp(x)
 * Allowed input range: anything
 */
double Math_Exp(double x) {
	return Exp2(x * LOG2E);
}

/************
 * Math_Log *
 ************/

/* Calculates the 3rd/3rd degree rational function LOG2 2524 listed in the
 * book's appendix.
 *
 * Associated math function: log_2(x)
 * Allowed input range: [0.5, 1]
 * Precision: 8.32
 */
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

	double P, Q;
	int i;

	P = A_P[3];
	for (i = 2; i >= 0; i--) {
		P *= x;
		P += A_P[i];
	}

	Q = A_Q[3];
	for (i = 2; i >= 0; i--) {
		Q *= x;
		Q += A_Q[i];
	}

	return P / Q;
}

/* Reduces the range of log_2(x) by using the property that
 *   log_2(x) = (x's exponent part) + log_2(x's mantissa part)
 * So, by manipulating the bits of the double-precision floating point number
 * one can reduce the range of the logarithm function.
 *
 * Associated math function: log_2(x)
 * Allowed input range: anything
 */
double Log2(double x) {
	union { double d; cc_uint64 i; } doi;
	int integer_part;

	if (x <= 0.0)
		return DOUBLE_NAN;

	doi.d = x;
	integer_part = (doi.i >> 52);
	integer_part -= 1023;

	doi.i |= (((cc_uint64) 1023) << 52);
	doi.i &= ~(((cc_uint64) 1024) << 52);

	return integer_part + Log2Stage1(doi.d);
}

/* Uses the property that
 *   log_e(x) = log_2(x) * log_e(2).
 *
 * Associated math function: log_e(x)
 * Allowed input range: anything
 */
double Math_Log(double x) {
	return Log2(x) * LOGE2;
}

/*****************
 * main function *
 *****************/

int main(void) {
	double step = 20.0/1000000.0;
	double x;

	for (x = -10.0; x <= 10.0; x += step) {
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

