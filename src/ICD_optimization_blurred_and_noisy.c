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
    struct TIFF_img input_img, noisy_blurred_img, MAP_est_img;
    double **img,**h,**y_img, **MAP_est_x;
    double *x, *y, *e;
    double noisy_pixel, sigma_x_sq, sigma_x_sq_hat, sigma_w_sq, temp, v, theta1, theta2, xi_update, cost1, cost2, total_cost;
    int32_t i, j, k, K, N, row, column;
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

    h = (double **)get_img(input_img.width,input_img.height,sizeof(double));
    for(i = 0; i < 5; i++){
        for(j = 0; j < 5; j++){
            h[i][j] = (3 - abs(2 - i)) * (3 - abs(2 - j));
        }
    }


    /* Blur the image and add noise to the image */
    y_img = (double **)get_img(input_img.width,input_img.height,sizeof(double));

    /* Convolution with the filter*/
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            y_img[i][j] = 0;
            for(int32_t kernel_i = -2; kernel_i <= 2; kernel_i++){
                for(int32_t kernel_j = -2; kernel_j <= 2; kernel_j++){
                    int32_t circ_i = (input_img.height + i + kernel_i) % input_img.height;
                    int32_t circ_j = (input_img.width + j + kernel_j) % input_img.width;

                    y_img[i][j] += (h[2 + kernel_i][2 + kernel_j] * img[circ_i][circ_j]) / 81.0;
                }
            }
        }
    }





    /* set up structure for output achromatic image */
    /* to allocate a full color image use type 'c' */
    get_TIFF ( &noisy_blurred_img, input_img.height, input_img.width, 'g' );

    /* Set seed for random noise generator */
    srandom2(1);

    /* Adding noise to the image */
    for ( i = 0; i < input_img.height; i++ ){
        for ( j = 0; j < input_img.width; j++ ) {
            noisy_pixel = y_img[i][j] + 4 * normal(); /* Add noise N(0, 4^2) to image */;
            noisy_blurred_img.mono[i][j] = MAX(MIN((int32_t)noisy_pixel, 255), 0);
        }
    }

    /* open noisy image file */
    if ( ( fp = fopen ( "noisy_blurred_img.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file ncpe_img.tif\n");
        exit ( 1 );
    }

    /* write noisy image */
    if ( write_TIFF ( fp, &noisy_blurred_img ) ) {
        fprintf ( stderr, "error writing TIFF file %s\n", argv[2] );
        exit ( 1 );
    }

    /* close noisy image file */
    fclose ( fp );

    /* ICD optimization */
    /* Set K = desired number of iterations */
    K = 20;

    /* Get the estimate sigma_x_sq_hat */
    sigma_x_sq_hat = 0;
    double p = 2;
    N = input_img.height * input_img.width;
    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            for(k = 0; k < 8; k++){
                temp = fabs(img[i][j] - img[(input_img.height + neighborhood[k].x + i) % input_img.height][(input_img.width + neighborhood[k].y + j) % input_img.width]);
                sigma_x_sq_hat += pow(temp, p) / neighborhood[k].g_inv;
            }
        }
    }

    sigma_x_sq_hat /= (double) N;
    
    sigma_x_sq_hat /= 2.0; // Here, we have to divide by 2 since we are considering all the unique unordered pairs (not ordered pairs)

    /* Select desired values of sigma_x and sigma_W */
    sigma_w_sq = 4.0 * 4.0;
    sigma_x_sq = sigma_x_sq_hat;

    /* Rasterize blurred and noise added image (y) and initialize x */
    int32_t len = input_img.height * input_img.width;
    x = (double *)calloc(len, sizeof(double));
    y = (double *)calloc(len, sizeof(double));

    for(i = 0; i < input_img.height; i++){
        for(j = 0; j < input_img.width; j++){
            y[i * input_img.width + j] = (double) noisy_blurred_img.mono[i][j];
            x[i * input_img.width + j] = (double) noisy_blurred_img.mono[i][j];
        }
    }


    /* Initial Cost Calculation */
    cost2 = 0.0;
    for(i = 0; i < len; i++){
        row = i / input_img.width;
        column = i % input_img.width;

        for(int32_t ind = 0; ind < 8; ind++){
            int32_t circ_i = (input_img.height + row + neighborhood[ind].x) % input_img.height;
            int32_t circ_j = (input_img.width + column + neighborhood[ind].y) % input_img.width;

            cost2 += pow(x[i] - x[circ_i * input_img.width + circ_j], 2) / neighborhood[ind].g_inv;
        }
    }

    cost2 /= 2 * 2 * sigma_x_sq;     // Here, we have to divide by an additional 2 since we are considering all the unique unordered pairs (not ordered pairs)

    cost1 = 0.0;
    /* Initialize e = y − Hx */
    e = (double *)calloc(len, sizeof(double));
    for(i = 0; i < len; i++){
        row = i / input_img.width;
        column = i % input_img.width;

        e[i] = y[i];
        for(int32_t kernel_i = -2; kernel_i <= 2; kernel_i++){
            for(int32_t kernel_j = -2; kernel_j <= 2; kernel_j++){
                int32_t circ_i = (input_img.height + row + kernel_i) % input_img.height;
                int32_t circ_j = (input_img.width + column + kernel_j) % input_img.width;

                e[i] -= (h[2 + kernel_i][2 + kernel_j] * x[circ_i * input_img.width + circ_j]) / 81.0;
            }
        }
        cost1 += e[i] * e[i];
    }

    cost1 /= 2 * sigma_w_sq;

    total_cost = cost1 + cost2;
    printf("%d %f\n", 0, total_cost);

    /* Perform iterations of ICD */
    for(k = 0; k < K; k++){
        cost1 = 0.0;
        cost2 = 0.0;
        for(i = 0; i < len; i++){
            row = i / input_img.width;
            column = i % input_img.width;

            v = x[i];    /* Update v */

            /* Update theta1 and theta2 */
            theta1 = 0.0;
            theta2 = 0.0;
            for(int32_t kernel_i = -2; kernel_i <= 2; kernel_i++){
                for(int32_t kernel_j = -2; kernel_j <= 2; kernel_j++){
                    // h is a space invariant filter
                    int32_t circ_i = (input_img.height + row + kernel_i) % input_img.height;
                    int32_t circ_j = (input_img.width + column + kernel_j) % input_img.width;

                    theta1 -= (e[circ_i * input_img.width + circ_j] * h[2 + kernel_i][2 + kernel_j]) / 81.0;
                    theta2 += (h[2 + kernel_i][2 + kernel_j] * h[2 + kernel_i][2 + kernel_j]) / (81.0 * 81.0);
                }
            }
            theta1 /= sigma_w_sq;
            theta2 /= sigma_w_sq;

            /* Update xi */
            xi_update = 0.0;
            for(int32_t ind = 0; ind < 8; ind++){
                int32_t circ_i = (input_img.height + row + neighborhood[ind].x) % input_img.height;
                int32_t circ_j = (input_img.width + column + neighborhood[ind].y) % input_img.width;

                xi_update += x[circ_i * input_img.width + circ_j] / neighborhood[ind].g_inv;
            }
            xi_update /= sigma_x_sq;
            xi_update = (xi_update + (theta2 * v) - theta1) / (theta2 + (1 / sigma_x_sq));
            x[i] = MAX(0.0, xi_update);

            /* Update e */
            for(int32_t kernel_i = -2; kernel_i <= 2; kernel_i++){
                for(int32_t kernel_j = -2; kernel_j <= 2; kernel_j++){
                    int32_t circ_i = (input_img.height + row + kernel_i) % input_img.height;
                    int32_t circ_j = (input_img.width + column + kernel_j) % input_img.width;

                    e[circ_i * input_img.width + circ_j] -= (h[2 + kernel_i][2 + kernel_j] * (x[i] - v)) / 81.0;
                }
            }
        }
        /* Cost Calculation */
        for(int i = 0; i < len; i++){
            cost1 += e[i] * e[i];

            row = i / input_img.width;
            column = i % input_img.width;

            for(int32_t ind = 0; ind < 8; ind++){
                int32_t circ_i = (input_img.height + row + neighborhood[ind].x) % input_img.height;
                int32_t circ_j = (input_img.width + column + neighborhood[ind].y) % input_img.width;

                cost2 += pow(x[i] - x[circ_i * input_img.width + circ_j], 2) / neighborhood[ind].g_inv;
            }
        }
        cost1 /= 2 * sigma_w_sq;
        cost2 /= 2 * 2 * sigma_x_sq;     // Here, we have to divide by an additional 2 since we are considering all the unique unordered pairs (not ordered pairs)
        total_cost = cost1 + cost2;
        printf("%d %f\n", k + 1, total_cost);
    }




    /* Allocate image of double precision floats */
    MAP_est_x = (double **)get_img(input_img.width,input_img.height,sizeof(double));

    /* Derasterize the array into an image */
    for(k = 0; k < len; k++){
        i = k / input_img.width;
        j = k % input_img.width;
        MAP_est_x[i][j] = x[k];
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
    if ( ( fp = fopen ( "MAP_est_blurred_noisy_img.tif", "wb" ) ) == NULL ) {
        fprintf ( stderr, "cannot open file MAP_est_blurred_noisy_img.tif\n");
        exit ( 1 );
    }

    /* write MAP estimate image */
    if ( write_TIFF ( fp, &MAP_est_img ) ) {
        fprintf ( stderr, "error writing TIFF file %s\n", argv[2] );
        exit ( 1 );
    }

    /* close MAP estimate image file */
    fclose ( fp );

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


