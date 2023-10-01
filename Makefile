
all: sin cos atan exp log libm_sin libm_cos libm_atan libm_exp libm_log

sin: approx.c
	gcc -O2 approx.c -DSIN -o sin

cos: approx.c
	gcc -O2 approx.c -DCOS -o cos

atan: approx.c
	gcc -O2 approx.c -DATAN -o atan

exp: approx.c
	gcc -O2 approx.c -DEXP -o exp

log: approx.c
	gcc -O2 approx.c -DLOG -o log

libm_sin: approx.c
	gcc -O2 approx.c -DLIBM_SIN -lm -o libm_sin

libm_cos: approx.c
	gcc -O2 approx.c -DLIBM_COS -lm -o libm_cos

libm_atan: approx.c
	gcc -O2 approx.c -DLIBM_ATAN -lm -o libm_atan

libm_exp: approx.c
	gcc -O2 approx.c -DLIBM_EXP -lm -o libm_exp

libm_log: approx.c
	gcc -O2 approx.c -DLIBM_LOG -lm -o libm_log

clean:
	rm -vf sin cos atan exp log libm_sin libm_cos libm_atan libm_exp libm_log sin_output.txt cos_output.txt atan_output.txt exp_output.txt log_output.txt libm_sin_output.txt libm_cos_output.txt libm_atan_output.txt libm_exp_output.txt libm_log_output.txt

test: all
	@echo "*** Running sin program"
	time ./sin > sin_output.txt
	@echo ""
	@echo "*** Running cos program"
	time ./cos > cos_output.txt
	@echo ""
	@echo "*** Running atan program"
	time ./atan > atan_output.txt
	@echo ""
	@echo "*** Running exp program"
	time ./exp > exp_output.txt
	@echo ""
	@echo "*** Running log program"
	time ./log > log_output.txt
	@echo ""
	@echo "*** Running libm_sin program"
	time ./libm_sin > libm_sin_output.txt
	@echo ""
	@echo "*** Running libm_cos program"
	time ./libm_cos > libm_cos_output.txt
	@echo ""
	@echo "*** Running libm_atan program"
	time ./libm_atan > libm_atan_output.txt
	@echo ""
	@echo "*** Running libm_exp program"
	time ./libm_exp > libm_exp_output.txt
	@echo ""
	@echo "*** Running libm_log program"
	time ./libm_log > libm_log_output.txt
	@echo ""
	@echo "*** Graphing"
	python3 graph.py

