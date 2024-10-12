#include <math.h>
#include <stdbool.h>
#include "tiff.h"
#include "allocate.h"
#include "randlib.h"
#include "typeutil.h"

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)

void error(char *name);

struct Neighbor{
    int32_t x;
    int32_t y;
    double g_inv;
};

int main (int argc, char **argv) {
    FILE *fp;
    struct TIFF_img input_img, noisy_img, MAP_est_img;
    double **img,**MAP_est_x;
    double noisy_pixel, sigma_x_sq, sigma_x_sq_hat, sigma_w_sq, temp, final_cost;
    int32_t i,j,k,d, K, N;
    struct Neighbor neighborhood[8];

    /* Define Neighborhood */
    k = 0;
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
    get_TIFF ( &noisy_img, input_img.height, input_img.width, 'g' );

    /* Set seed for random noise generator */
    srandom2(1);

    /* Adding noise to the image */
    for ( i = 0; i < input_img.height; i++ ){
        for ( j = 0; j < input_img.width; j++ ) {
            noisy_pixel = img[i][j] + 16.0 * normal(); /* Add noise N(0, 16^2) to image */;
            noisy_img.mono[i][j] = MAX(MIN((int32_t)noisy_pixel, 255), 0);
        }
    }

    /* open noisy image file */
    if ( ( fp = fopen ( "noisy_img.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file ncpe_img.tif\n");
        exit ( 1 );
    }

    /* write noisy image */
    if ( write_TIFF ( fp, &noisy_img ) ) {
        fprintf ( stderr, "error writing TIFF file %s\n", argv[2] );
        exit ( 1 );
    }

    /* close noisy image file */
    fclose ( fp );

    K = 20;

    /* Get the estimate sigma_x_sq_hat */
    sigma_x_sq_hat = 0;
    N = input_img.height * input_img.width;
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            for(k = 0; k < 8; k++){
                temp = fabs(img[i][j] - img[(input_img.height + neighborhood[k].x + i) % input_img.height][(input_img.width + neighborhood[k].y + j) % input_img.width]);
                sigma_x_sq_hat += pow(temp, 2) / neighborhood[k].g_inv;
            }
        }
    }

    sigma_x_sq_hat /= (double) N;

    sigma_w_sq = 16.0 * 16.0;
    sigma_x_sq = sigma_x_sq_hat;

    printf("%f %f\n",sigma_w_sq, sigma_x_sq);

    /* Allocate image of double precision floats */
    MAP_est_x = (double **)get_img(input_img.width,input_img.height,sizeof(double));

    /* For each i ∈ S, initialize with ML estimate */
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            MAP_est_x[i][j] = noisy_img.mono[i][j];
        }
    }

    double cost1 = 0.0;
    double cost2 = 0.0;
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            for(k = 0; k < 8; k++){
                cost2 += pow(MAP_est_x[i][j] - MAP_est_x[(input_img.height + neighborhood[k].x + i) % input_img.height][(input_img.width + neighborhood[k].y + j) % input_img.width], 2) / neighborhood[k].g_inv;
            }
        }
    }

    final_cost = cost2 / (2.0 * sigma_x_sq);

    FILE *file = fopen("costs_sigma_1.txt", "w");
    if (file == NULL) {
        perror("Error opening file");
        return -1;
    }
    fprintf(file, "%d %f\n", 0, final_cost);

    for (k = 0; k < K; k++){
        cost1 = 0.0;
        cost2 = 0.0;
        for(i = 0; i < input_img.height; i++){
            for(j = 0; j < input_img.width; j++){
                temp = 0;
                for(d = 0; d < 8; d++){
                    double x_j = MAP_est_x[(input_img.height + neighborhood[d].x + i) % input_img.height][(input_img.width + neighborhood[d].y + j) % input_img.width];
                    temp += x_j / neighborhood[d].g_inv;
                }
                temp *= (sigma_w_sq / sigma_x_sq);
                temp += (double) noisy_img.mono[i][j];
                temp /= (1 + (sigma_w_sq / sigma_x_sq));
                // printf("%f %f\n", MAP_est_x[i][j], temp);
                MAP_est_x[i][j] = MAX(0, temp);

                /* Calculate cost 1*/
                cost1 += pow((double) noisy_img.mono[i][j] - MAP_est_x[i][j], 2);

                /* Calculate cost 2*/
                for(d = 0; d < 8; d++){
                    cost2 += pow(MAP_est_x[i][j] - MAP_est_x[(input_img.height + neighborhood[d].x + i) % input_img.height][(input_img.width + neighborhood[d].y + j) % input_img.width], 2) / neighborhood[d].g_inv;
                }
            }
        }
        final_cost = (cost1 / (2 * sigma_w_sq)) + (cost2 / (2 * sigma_x_sq));
        fprintf(file, "%d %f\n", k + 1, final_cost);
    }

    /* set up structure for output achromatic image */
    /* to allocate a full color image use type 'c' */
    get_TIFF ( &MAP_est_img, input_img.height, input_img.width, 'g' );

    /* Print out the resulting MAP estimate */
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            MAP_est_img.mono[i][j] = (int32_t)MIN(255.0, MAP_est_x[i][j]);
        }
    }

    /* open MAP estimate image file */
    if ( ( fp = fopen ( "MAP_est_img_sigma_x_1.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file MAP_est_img.tif\n");
        exit ( 1 );
    }

    /* write MAP estimate image */
    if ( write_TIFF ( fp, &MAP_est_img ) ) {
        fprintf ( stderr, "error writing TIFF file %s\n", argv[2] );
        exit ( 1 );
    }

    /* close MAP estimate image file */
    fclose ( fp );


    /******************************* Let's change the value of sigma_x *******************************/
    sigma_x_sq = 5 * sigma_x_sq_hat;

    printf("%f %f\n",sigma_w_sq, sigma_x_sq);

    /* For each i ∈ S, initialize with ML estimate */
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            MAP_est_x[i][j] = noisy_img.mono[i][j];
        }
    }

    /* Calculate the initial cost */
    cost1 = 0.0;
    cost2 = 0.0;
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            for(k = 0; k < 8; k++){
                cost2 += pow(MAP_est_x[i][j] - MAP_est_x[(input_img.height + neighborhood[k].x + i) % input_img.height][(input_img.width + neighborhood[k].y + j) % input_img.width], 2) / neighborhood[k].g_inv;
            }
        }
    }

    final_cost = cost2 / (2.0 * sigma_x_sq);

    file = fopen("costs_sigma_5.txt", "w");
    if (file == NULL) {
        perror("Error opening file");
        return -1;
    }
    fprintf(file, "%d %f\n", 0, final_cost);

    for (k = 0; k < K; k++){
        cost1 = 0.0;
        cost2 = 0.0;
        for(i = 0; i < input_img.height; i++){
            for(j = 0; j < input_img.width; j++){
                temp = 0;
                for(d = 0; d < 8; d++){
                    double x_j = MAP_est_x[(input_img.height + neighborhood[d].x + i) % input_img.height][(input_img.width + neighborhood[d].y + j) % input_img.width];
                    temp += x_j / neighborhood[d].g_inv;
                }
                temp *= (sigma_w_sq / sigma_x_sq);
                temp += (double) noisy_img.mono[i][j];
                temp /= (1 + (sigma_w_sq / sigma_x_sq));

                MAP_est_x[i][j] = MAX(0, temp);

                /* Calculate cost 1*/
                cost1 += pow((double) noisy_img.mono[i][j] - MAP_est_x[i][j], 2);

                /* Calculate cost 2*/
                for(d = 0; d < 8; d++){
                    cost2 += pow(MAP_est_x[i][j] - MAP_est_x[(input_img.height + neighborhood[d].x + i) % input_img.height][(input_img.width + neighborhood[d].y + j) % input_img.width], 2) / neighborhood[d].g_inv;
                }
            }
        }
        final_cost = (cost1 / (2 * sigma_w_sq)) + (cost2 / (2 * sigma_x_sq));
        fprintf(file, "%d %f\n", k + 1, final_cost);
    }


    /* set up structure for output achromatic image */
    /* to allocate a full color image use type 'c' */
    get_TIFF ( &MAP_est_img, input_img.height, input_img.width, 'g' );

    /* Print out the resulting MAP estimate */
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            MAP_est_img.mono[i][j] = (int32_t)MIN(255.0, MAP_est_x[i][j]);
        }
    }

    /* open MAP estimate image file */
    if ( ( fp = fopen ( "MAP_est_img_sigma_x_5.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file MAP_est_img_sigma_x_5.tif\n");
        exit ( 1 );
    }

    /* write MAP estimate image */
    if ( write_TIFF ( fp, &MAP_est_img ) ) {
        fprintf ( stderr, "error writing TIFF file %s\n", argv[2] );
        exit ( 1 );
    }

    /* close MAP estimate image file */
    fclose ( fp );


    /******************************* Let's change the value of sigma_x_sq -> sigma_x_sq = (1/5) * sigma_x_sq_hat *******************************/
    sigma_x_sq = sigma_x_sq_hat / 5;

    printf("%f %f\n",sigma_w_sq, sigma_x_sq);

    /* For each i ∈ S, initialize with ML estimate */
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            MAP_est_x[i][j] = noisy_img.mono[i][j];
        }
    }

    /* Calculate the initial cost */
    cost1 = 0.0;
    cost2 = 0.0;
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            for(k = 0; k < 8; k++){
                cost2 += pow(MAP_est_x[i][j] - MAP_est_x[(input_img.height + neighborhood[k].x + i) % input_img.height][(input_img.width + neighborhood[k].y + j) % input_img.width], 2) / neighborhood[k].g_inv;
            }
        }
    }

    final_cost = cost2 / (2.0 * sigma_x_sq);

    file = fopen("costs_sigma_1_5.txt", "w");
    if (file == NULL) {
        perror("Error opening file");
        return -1;
    }
    fprintf(file, "%d %f\n", 0, final_cost);

    for (k = 0; k < K; k++){
        cost1 = 0.0;
        cost2 = 0.0;
        for(i = 0; i < input_img.height; i++){
            for(j = 0; j < input_img.width; j++){
                temp = 0;
                for(d = 0; d < 8; d++){
                    double x_j = MAP_est_x[(input_img.height + neighborhood[d].x + i) % input_img.height][(input_img.width + neighborhood[d].y + j) % input_img.width];
                    temp += x_j / neighborhood[d].g_inv;
                }
                temp *= (sigma_w_sq / sigma_x_sq);
                temp += (double) noisy_img.mono[i][j];
                temp /= (1 + (sigma_w_sq / sigma_x_sq));

                MAP_est_x[i][j] = MAX(0, temp);

                /* Calculate cost 1*/
                cost1 += pow((double) noisy_img.mono[i][j] - MAP_est_x[i][j], 2);

                /* Calculate cost 2*/
                for(d = 0; d < 8; d++){
                    cost2 += pow(MAP_est_x[i][j] - MAP_est_x[(input_img.height + neighborhood[d].x + i) % input_img.height][(input_img.width + neighborhood[d].y + j) % input_img.width], 2) / neighborhood[d].g_inv;
                }
            }
        }
        final_cost = (cost1 / (2 * sigma_w_sq)) + (cost2 / (2 * sigma_x_sq));
        fprintf(file, "%d %f\n", k + 1, final_cost);
    }


    /* set up structure for output achromatic image */
    /* to allocate a full color image use type 'c' */
    get_TIFF ( &MAP_est_img, input_img.height, input_img.width, 'g' );

    /* Print out the resulting MAP estimate */
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            MAP_est_img.mono[i][j] = (int32_t)MIN(255.0, MAP_est_x[i][j]);
        }
    }

    /* open MAP estimate image file */
    if ( ( fp = fopen ( "MAP_est_img_sigma_x_1_5.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file MAP_est_img_sigma_x_1_5.tif\n");
        exit ( 1 );
    }

    /* write MAP estimate image */
    if ( write_TIFF ( fp, &MAP_est_img ) ) {
        fprintf ( stderr, "error writing TIFF file %s\n", argv[2] );
        exit ( 1 );
    }

    /* close MAP estimate image file */
    fclose ( fp );


    // standard_ICD(K, &input_img, &noisy_img, cost1, cost2, sigma_x_sq, sigma_w_sq, file, MAP_est_x, neighborhood);
    

    // printf("Hello Praditha!");
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


