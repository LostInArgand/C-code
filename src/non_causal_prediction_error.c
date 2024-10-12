#include <math.h>
#include <stdbool.h>
#include "tiff.h"
#include "allocate.h"
#include "randlib.h"
#include "typeutil.h"

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)

void error(char *name);
bool valid(int32_t x, int32_t y, int32_t N, int32_t M);

struct Neighbor{
    int32_t x;
    int32_t y;
    double g_inv;
};

int main (int argc, char **argv) {
    FILE *fp;
    struct TIFF_img input_img, green_img, color_img, ncpe_img;
    double **img;
    double ncpe_pixel;
    int32_t i,j,pixel;
    struct Neighbor neighborhood[8];

    /* Define Neighborhood */
    int k = 0;
    for(i=-1; i<=1; i++){
        for(j=-1; j<=1; j++){
            if(i==0 && j==0) continue; // Pixel itself is not a neighbor
            neighborhood[k].x = i;
            neighborhood[k].y = j;
            neighborhood[k].g_inv = 6.0 * (abs(i) + abs(j));
            k++;
        }
    }

    if ( argc != 2 ) error( argv[0] );

    /* open image file */
    if ( ( fp = fopen ( argv[1], "rb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file %s\n", argv[1] );
        exit ( 1 );
    }

    /* read image */
    if ( read_TIFF ( fp, &input_img ) ) {
        fprintf ( stderr, "error reading file %s\n", argv[1] );
        exit ( 1 );
    }

    /* close image file */
    fclose ( fp );

    /* check the type of image data */
    if ( input_img.TIFF_type != 'g' ) {
        fprintf ( stderr, "error:  image must be grayscaled image\n" );
        exit ( 1 );
    }

    /* Allocate image of double precision floats */
    img = (double **)get_img(input_img.width,input_img.height,sizeof(double));

    /* copy grayscaled component to double array */
    for ( i = 0; i < input_img.height; i++ ){
        for ( j = 0; j < input_img.width; j++ ) {
            img[i][j] = input_img.mono[i][j];
        }
    }


    /* set up structure for output achromatic image */
    /* to allocate a full color image use type 'c' */
    get_TIFF ( &ncpe_img, input_img.height, input_img.width, 'g' );

    /* Computing the non causal prediction filter */
    for ( i = 0; i < input_img.height; i++ ){
        for ( j = 0; j < input_img.width; j++ ) {
            ncpe_pixel = img[i][j];
            for(int32_t k = 0; k < 8; k++){
                ncpe_pixel -= img[(input_img.height + neighborhood[k].x + i) % input_img.height][(input_img.width + neighborhood[k].y + j) % input_img.width] / neighborhood[k].g_inv;
            }
            pixel = (int32_t)(ncpe_pixel + 127.0);
            ncpe_img.mono[i][j] = MAX(MIN(pixel, 255), 0);
        }
    }

    /* open green image file */
    if ( ( fp = fopen ( "ncpe_img.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file ncpe_img.tif\n");
        exit ( 1 );
    }

    /* write green image */
    if ( write_TIFF ( fp, &ncpe_img ) ) {
        fprintf ( stderr, "error writing TIFF file %s\n", argv[2] );
        exit ( 1 );
    }

    /* close green image file */
    fclose ( fp );

    /* ML estimator for p values between 0.1 to 2*/
    double p = 0.1, dp = 0.01, sigma;

    /* allocate memory */
    int32_t ind = 0;
    int32_t N = input_img.height * input_img.width;

    while(p < 2.001){
        sigma = 0;
        for ( i = 0; i < input_img.height; i++ ){
            for ( j = 0; j < input_img.width; j++ ){
                for(int32_t k = 0; k < 8; k++){
                    /* avoid boundary conditions */
                    if(valid(neighborhood[k].x + i, neighborhood[k].y + j, input_img.height, input_img.width)){
                        pixel = fabs(img[i][j] - img[neighborhood[k].x + i][neighborhood[k].y + j]);
                        pixel = pow(pixel, p);
                        sigma += pixel / neighborhood[k].g_inv;
                        // fprintf(stdout, "%g\n", sigma);
                    }
                }
            }
        }
        sigma /= N;
        sigma = pow(sigma, 1.0/p);

        /* Output the sigma and p values*/
        fprintf ( stdout, "%g, %g; \n", p, sigma );
        p += dp;
    }

    return(0);
}

void error(char *name)
{
    printf("usage:  %s  image.tiff \n\n",name);
    printf("this program reads in a 24-bit color TIFF image.\n");
    printf("It then horizontally filters the green component, adds noise,\n");
    printf("and writes out the result as an 8-bit image\n");
    printf("with the name 'green.tiff'.\n");
    printf("It also generates an 8-bit color image,\n");
    printf("that swaps red and green components from the input image");
    exit(1);
}

bool valid(int32_t x, int32_t y, int32_t N, int32_t M){
    if(0 > x || x >= N) return false;
    if(0 > y || y >= M) return false;
    return true;
}

