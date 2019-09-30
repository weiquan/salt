all : snpaln snpindex polish
snpaln:
	cd ./Align_src/ && $(MAKE)
	mv ./Align_src/snpaln ./Bin/
	mv ./Align_src/snpaln_debug ./Debug/

snpindex:
	cd ./Index_src/ && $(MAKE)
	mv ./Index_src/index ./Bin/
	mv ./Index_src/index_debug ./Debug/
polish:
	cd ./Polish_src/ && $(MAKE)
	mv ./Polish_src/polish ./Bin/
clean:
	cd ./Align_src/ && make clean
	cd ./Index_src/ && make clean
	cd ./Polish_src/ && make clean
	
	rm -f ./Bin/*
	rm -f ./Debug/*
	#rm -f ./Test/Index/*
	#rm -f ./Test/Reads/*
	#rm -f ./Test/SAM/*
	#rm -f ./Test/Variants/*
