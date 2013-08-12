/*
 * Plot.h
 *
 *  Created on: Jun 26, 2013
 *      Author: Marc Neveu (mneveu@asu.edu)
 *
 *  Useful SDL 2.0 functions from:
 *  http://twinklebeardev.blogspot.com/2012/07/lesson-2-dont-put-everything-in-main.html
 */

#ifndef PLOT_H_
#define PLOT_H_

#include <stdio.h>
#include <string.h>

#include <SDL.h>
#include <SDL_image.h>
#include <SDL_mixer.h>

const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 600;

SDL_Texture* LoadImage(char * filename);
void ApplySurface(int x, int y, SDL_Texture *tex, SDL_Renderer *rend, SDL_Rect *clip);

SDL_Window *window = NULL;
SDL_Renderer *renderer = NULL;

// A function for loading textures from a filename.
// We can pass a file name as a string and get back a pointer to the loaded SDL_Texture.
// Note that the pointer will be nullptr if the loading fails because we initialize both pointers as NULL for error checking.

SDL_Texture* LoadImage(char * filename) {
	SDL_Texture* tex = NULL;
	tex = IMG_LoadTexture(renderer, filename);
	if (tex == NULL)
		printf("Plot: Error: Could not load image\n");
    return tex;
}

// A function to simplify our draw calls and also allow us to specify a position to draw the image too on the screen.
// We'll want it to be able to take an x, y coordinate position along with a texture pointer and a renderer pointer
// and then draw the texture to that position.

void ApplySurface(int x, int y, SDL_Texture *tex, SDL_Renderer *rend, SDL_Rect *clip) {
    SDL_Rect pos;
    pos.x = x;
    pos.y = y;
    //Detect if we should use clip width settings or texture width
    if (clip != NULL){
        pos.w = clip->w;
        pos.h = clip->h;
    }
    else {
        SDL_QueryTexture(tex, NULL, NULL, &pos.w, &pos.h);
    }
    SDL_RenderCopy(rend, tex, clip, &pos);
}

// Clip numbers from the Numbers.png texture.

SDL_Rect ClipNumber(int number,int font) {
	SDL_Rect clip;
	if (font == 12) {
		clip.h = 20;
		clip.w = 12;
		clip.y = 4;
		clip.x = number*10 + 4;
	}
	else if (font == 14) {
		clip.h = 20;
		clip.w = 12;
		clip.y = 22;
		clip.x = number*10 + 4;
	}
	else if (font == 18) {
		clip.h = 20;
		if (number == 1) clip.w = 10;
		else clip.w = 12;
		clip.y = 46;
		clip.x = number*10 + 4;
	}
	else {
		printf("Plot: ClipNumber: Font size not supported.\n");
		clip.x = 0, clip.y = 0, clip.w = 0, clip.h = 0;
	}
	return clip;
}

// Get a pixel

Uint32 get_pixel32( SDL_Surface *surface, int x, int y ) {
    //Convert the pixels to 32 bit
    Uint32 *pixels = (Uint32 *)surface->pixels;

    //Get the requested pixel
    return pixels[ ( y * surface->w ) + x ];
}

// Put a pixel

void put_pixel32( SDL_Surface *surface, int x, int y, Uint32 pixel ) {
    //Convert the pixels to 32 bit
    Uint32 *pixels = (Uint32 *)surface->pixels;

    //Set the pixel
    pixels[ ( y * surface->w ) + x ] = pixel;
}

#endif /* PLOT_H_ */
