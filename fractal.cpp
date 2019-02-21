// Fractal.cpp : Defines the entry point for the console application.
//

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <complex>
#include <iostream>
#include <time.h>
#include <mpi.h>

#define NROWS_TAG 0
#define DATA_TAG 1

void WriteTGA_RGB(const char* filename, unsigned char* data, unsigned int width, unsigned int height)
{
	FILE *f = fopen(filename, "wb");
	if (!f) {
		fprintf(stderr, "Unable to create output TGA image `%s'\n", filename);
		exit(EXIT_FAILURE);
	}

	fputc(0x00, f); /* ID Length, 0 => No ID        */
	fputc(0x00, f); /* Color Map Type, 0 => No color map included   */
	fputc(0x02, f); /* Image Type, 2 => Uncompressed, True-color Image */
	fputc(0x00, f); /* Next five bytes are about the color map entries */
	fputc(0x00, f); /* 2 bytes Index, 2 bytes length, 1 byte size */
	fputc(0x00, f);
	fputc(0x00, f);
	fputc(0x00, f);
	fputc(0x00, f); /* X-origin of Image    */
	fputc(0x00, f);
	fputc(0x00, f); /* Y-origin of Image    */
	fputc(0x00, f);
	fputc(width & 0xff, f); /* Image Width      */
	fputc((width >> 8) & 0xff, f);
	fputc(height & 0xff, f); /* Image Height     */
	fputc((height >> 8) & 0xff, f);
	fputc(0x18, f); /* Pixel Depth, 0x18 => 24 Bits */
	fputc(0x20, f); /* Image Descriptor     */

	for (int y = height - 1; y >= 0; y--) {
		for (size_t x = 0; x < width; x++) {
			const size_t i = (y * width + x) * 3;
			fputc(data[i + 2], f); /* write blue */
			fputc(data[i + 1], f); /* write green */
			fputc(data[i], f); /* write red */
		}
	}
}

int main(int argc, char **argv)
{

	const unsigned int domainWidth = 1024;
	const unsigned int domainHeight = 1024;
	unsigned char *data = new unsigned char[domainWidth * domainHeight * 3];
	std::memset(data, 0, domainWidth * domainHeight * 3 * sizeof(unsigned char));

	std::complex<double> c(0.285, 0.013);
	std::complex<double> K(0.353, 0.288);
	std::complex<double> center(-1.23, -1.23);
	double scale = 2.35;

	double startTime, endTime;

	std::cout.precision(12);

	const unsigned int maxIterations = 100;


	//variables for MPI communication:
	int id, nproc;
	MPI_Status status;

	// Initialize MPI:
	MPI_Init(&argc, &argv);
	// Get my rank:
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	// Get the total number of processors:
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
 	
	MPI_Barrier(MPI_COMM_WORLD) ;//for precize timing


	unsigned int nRows;
	unsigned int temp_nRows;
	unsigned int startingRow;
	unsigned int split = domainHeight / (nproc-1);
	unsigned int extraLines = domainHeight % (nproc-1);
	unsigned int position;
	unsigned int temp_position;

	
	if(id!=0){

		if(id <= extraLines) {
			nRows = split + 1;
			startingRow = (id-1) * (split+1);
		}else{
			nRows = split;
			startingRow = (extraLines * (split+1)) + (((id-1) - extraLines) * split);
		}


		MPI_Send(&nRows, 1, MPI_INT, 0, NROWS_TAG, MPI_COMM_WORLD);
		std::cout<<"I am the node "<<id<<"; rows: "<<nRows<<std::endl;

		unsigned int pSize = nRows*domainWidth*maxIterations*3;
		unsigned int pixels[pSize];

		for (unsigned int y = startingRow; y < startingRow + nRows; y++) {
			for (unsigned int x = 0; x < domainWidth; x++) {
				std::complex<double> z(x / (double)domainWidth * scale + center.real(),
					y / (double)domainHeight * scale + center.imag());

				for (unsigned int iteration = 0; iteration < maxIterations; ++iteration)
				{
					z = z * z + K;

					
					if (std::abs(z) > 1.0f)
					{
						//send here the coordinate
						//pixels[(x + y * domainWidth) * 3 +0] = 255;
						//pixels[(x + y * domainWidth) * 3 +1] = 255;
						//pixels[(x + y * domainWidth) * 3 +2] = 255;
						
						
					}

				}
			}
			
		}
		std::cout<<"I am the node "<<id<<"sizeof(pixels)"<<sizeof(pixels)<<std::endl;
		MPI_Send(&pixels, pSize, MPI_INT, 0, DATA_TAG, MPI_COMM_WORLD);
		std::cout<<"1"<<std::endl;
	}
	if(id == 0){ //Master
		std::cout<<"2"<<std::endl;
		for(int i_id = 1; i_id < nproc; i_id++){

			MPI_Recv(&temp_nRows, 1, MPI_INT, i_id, NROWS_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			std::cout<<"temp_nRows:"<<temp_nRows<<std::endl;
			
			int pSize = temp_nRows*domainWidth*maxIterations*3;
			unsigned int pixels[pSize];
			std::cout<<"I am the node "<<i_id<<"sizeof(pixels)"<<sizeof(pixels)<<std::endl;
			MPI_Recv(&pixels, pSize, MPI_INT, i_id, DATA_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			std::cout<<"I am the node "<<i_id<<"sizeof(pixels)"<<sizeof(pixels)<<std::endl;
			
			



      		
		}
		std::cout<<"5"<<std::endl;
		//WriteTGA_RGB("julia.tga", data, domainWidth, domainHeight);
	}
	//startTime=MPI_Wtime(); // we get a time value at this point

//	for (unsigned int y = 0; y < domainHeight; ++y)
//	{
//		for (unsigned int x = 0; x < domainWidth; ++x)
//		{
//			std::complex<double> z(x / (double)domainWidth * scale + center.real(),
//				y / (double)domainHeight * scale + center.imag());
//
//			//std::complex<double> z(c);
//
//			for (unsigned int iteration = 0; iteration < maxIterations; ++iteration)
//			{
//				z = z * z + K;
//				if (std::abs(z) > 1.0f)
//				{
//					data[(x + y * domainWidth) * 3 + 0] = 255;
//					data[(x + y * domainWidth) * 3 + 1] = 255;
//					data[(x + y * domainWidth) * 3 + 2] = 255;
//				}
//			}
//		}
//	}
//
//	endTime = omp_get_wtime();
//
//	std::cout << "Computation time: " << endTime - startTime << std::endl;
//	std::cout << "Precision: " << omp_get_wtick() << std::endl;

	/* All processes clean up the MPI environment */
  	MPI_Finalize();
	delete[] data;
	return 0;
}
