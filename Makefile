all: 
#	@cd Lz4 ; make lz4
	@cd Hipo; make
	@cd Banks; make
#	@cd Utils; make

clean:
#	@cd Lz4; make clean
	@cd Hipo; make clean
	@cd Banks; make clean
#	@cd Utils; make clean
	@echo "Cleaning lib directory"
	@rm -rf lib/*
