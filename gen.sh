#!/bin/bash
mkdir -p release &&
	cd release &&
	cmake -DCMAKE_BUILD_TYPE=Release .. &&
	make &&
	test/intervaltest &&
	ln -sv ../etc/datasets.csv . &&
../etc/benchmark.sh


