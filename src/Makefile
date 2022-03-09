include ../make.cfg/flags.mak

all :	
	g++ pacer_imager.cpp -c $(OPT) -I $(BIGHORNS)/software/common -I $(BIGHORNS)/software/analysis/fitslib/ -D_UNIX
	g++ antenna_positions.cpp -c $(OPT) -I $(BIGHORNS)/software/common -I../fitslib/ -D_UNIX
	g++ pacer_imager_main.cpp pacer_imager.o antenna_positions.o -o pacer_dirty_image $(OPT) -I $(BIGHORNS)/software/common -I $(BIGHORNS)/software/analysis/fitslib/ $(COMMON_LIBS) -D_UNIX -lfftw3	
	cp pacer_dirty_image $(BIGHORNS)/bin/
	

test : old
	
old :
#	g++ pacer_imager.cpp -c $(OPT) -I $(BIGHORNS)/software/common -I../fitslib/ -D_UNIX
	g++ main.cpp -o pacer_dirty_image_old $(OPT) -I $(BIGHORNS)/software/common -I../fitslib/ $(COMMON_LIBS) -D_UNIX -lfftw3
	g++ main_test.cpp -o pacer_dirty_image_test $(OPT) -I $(BIGHORNS)/software/common -I../fitslib/ $(COMMON_LIBS) -D_UNIX -lfftw3
	g++ generator.cpp -o pacer_uvgrid_generator $(OPT) -I $(BIGHORNS)/software/common -I../fitslib/ $(COMMON_LIBS) -D_UNIX -lfftw3
	
	cp pacer_dirty_image_old $(BIGHORNS)/bin/
	cp pacer_dirty_image_test $(BIGHORNS)/bin/
	cp pacer_uvgrid_generator $(BIGHORNS)/bin/
	
clean :
	rm -f pacer_dirty_image pacer_uvgrid_generator *.o
		