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
#include <SDL_ttf.h>

const int SCREEN_WIDTH = 800;
const int SCREEN_HEIGHT = 600;

SDL_Texture* LoadImage(char * filename);
int ApplySurface(int x, int y, SDL_Texture *tex, SDL_Renderer *rend, SDL_Rect *clip);
SDL_Rect ClipNumber(int number,int font);
Uint32 get_pixel32( SDL_Surface *surface, int x, int y );               // Currently not used
int put_pixel32( SDL_Surface *surface, int x, int y, Uint32 pixel );   // Currently not used
int renderTexture(SDL_Texture *tex, SDL_Renderer *ren, int x, int y);
int scanNumber(char (*nb)[20], double i);
int File2tex(char* title, SDL_Texture** tex, char* path);
int File2surf(char* title, SDL_Surface** surf, char* path);
int darkenscreen(SDL_Surface **surf, SDL_Texture *bgtex, char FontFile[1024], char text[40], SDL_Renderer *renderer, int i,	int imax);

SDL_Window* window;
SDL_Renderer* renderer;

//-------------------------------------------------------------------
//                          Load Image
//-------------------------------------------------------------------

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

//-------------------------------------------------------------------
//                          Apply surface
//-------------------------------------------------------------------

// A function to simplify our draw calls and also allow us to specify a position to draw the image too on the screen.
// We'll want it to be able to take an x, y coordinate position along with a texture pointer and a renderer pointer
// and then draw the texture to that position.

int ApplySurface(int x, int y, SDL_Texture *tex, SDL_Renderer *rend, SDL_Rect *clip) {
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

    return 0;
}

//-------------------------------------------------------------------
//                           Get a pixel
//-------------------------------------------------------------------

Uint32 get_pixel32( SDL_Surface *surface, int x, int y ) {
    //Convert the pixels to 32 bit
    Uint32 *pixels = (Uint32 *)surface->pixels;

    //Get the requested pixel
    return pixels[ ( y * surface->w ) + x ];
}

//-------------------------------------------------------------------
//                           Put a pixel
//-------------------------------------------------------------------

int put_pixel32( SDL_Surface *surface, int x, int y, Uint32 pixel ) {
    //Convert the pixels to 32 bit
    Uint32 *pixels = (Uint32 *)surface->pixels;

    //Set the pixel
    pixels[ ( y * surface->w ) + x ] = pixel;

    return 0;
}

/**
* Render the message we want to display to a texture for drawing
* @param message The message we want to display
* @param fontFile The font we want to use to render the text
* @param color The color we want the text to be
* @param fontSize The size we want the font to be
* @param renderer The renderer to load the texture in
* @return An SDL_Texture containing the rendered message, or NULL if something went wrong
*/
SDL_Texture* renderText(char* message, char* fontFile, SDL_Color color, int fontSize, SDL_Renderer *renderer) {

	// Open the font
	TTF_Font *font = TTF_OpenFont(fontFile, fontSize);
	if (font == NULL){
		printf("IcyDwarf: Plot: renderText: Font file not loaded.\n");
		return NULL;
	}
	// First render to a surface as that's what TTF_RenderText returns, then load that surface into a texture
	SDL_Surface *surf = TTF_RenderText_Blended(font, message, color);
	if (surf == NULL){
		TTF_CloseFont(font);
		printf("IcyDwarf: Plot: renderText: Font not rendered to surface.\n");
		return NULL;
	}
	SDL_Texture *texture = SDL_CreateTextureFromSurface(renderer, surf);
	if (texture == NULL){
		printf("IcyDwarf: Plot: renderText: Font surface not converted to texture.\n");
	}
	// Clean up the surface and font
	SDL_FreeSurface(surf);
	TTF_CloseFont(font);
	return texture;
}

/**
* Draw an SDL_Texture to an SDL_Renderer at position x, y, preserving
* the texture's width and height
* @param tex The source texture we want to draw
* @param ren The renderer we want to draw too
* @param x The x coordinate to draw too
* @param y The y coordinate to draw too
*/
int renderTexture(SDL_Texture *tex, SDL_Renderer *ren, int x, int y) {
	//Setup the destination rectangle to be at the position we want
	SDL_Rect dst;
	dst.x = x;
	dst.y = y;
	//Query the texture to get its width and height to use
	SDL_QueryTexture(tex, NULL, NULL, &dst.w, &dst.h);
	SDL_RenderCopy(ren, tex, NULL, &dst);

	return 0;
}

//-------------------------------------------------------------------
//             Convert number to string, right-aligned
//-------------------------------------------------------------------

int scanNumber(char (*nb)[20], double i) {
	if (i<100.0) sprintf((*nb), "   %g", i);
	else if (i<1000.0) sprintf((*nb), "  %g", i);
	else sprintf((*nb), "%g", i);
	return 0;
}

//-------------------------------------------------------------------
//                 Open image file into a texture
//-------------------------------------------------------------------

int File2tex(char* title, SDL_Texture** tex, char* path) {
	char *temp_png = (char*)malloc(1024); // Don't forget to free!
	temp_png[0] = '\0';
	if (v_release == 1) strncat(temp_png,path,strlen(path)-20);
	else if (cmdline == 1) strncat(temp_png,path,strlen(path)-22);
	strcat(temp_png,title);
	(*tex) = LoadImage(temp_png);
	if ((*tex) == NULL) printf("Plot: Layer %s not loaded:\n",title);
	free(temp_png);

	return 0;
}

//-------------------------------------------------------------------
//                 Open image file into a surface
//-------------------------------------------------------------------

int File2surf(char* title, SDL_Surface** surf, char* path) {
	char *temp_png = (char*)malloc(1024); // Don't forget to free!
	temp_png[0] = '\0';
	if (v_release == 1) strncat(temp_png,path,strlen(path)-20);
	else if (cmdline == 1) strncat(temp_png,path,strlen(path)-22);
	strcat(temp_png,title);
	(*surf) = IMG_Load(temp_png);
	if ((*surf) == NULL) printf("Plot: Layer %s not loaded:\n",title);
	free(temp_png);

	return 0;
}

//-------------------------------------------------------------------
//                         Darken screen
//-------------------------------------------------------------------

int darkenscreen(SDL_Surface **surf, SDL_Texture *bgtex, char FontFile[1024], char text[40], SDL_Renderer *renderer, int i,	int imax) {

	SDL_Texture *loading_tex= NULL;
	SDL_Texture *surf_tex = NULL;
	SDL_Color white;
	white.r = 255; white.g = 255; white.b = 255; white.a = 50;

	Uint32 *pixmem32;
	int x = 0; int y = 0;

	for (x=0;x<(*surf)->w;x++) {
		for (y=0;y<=(*surf)->h;y++) {
			if ((y==290 || y==295) && (x>350 && x<=450)) {
				pixmem32 = (Uint32*) (*surf)->pixels + y*(*surf)->w + x;
				*pixmem32 = SDL_MapRGBA((*surf)->format, 255, 255, 255, 0);
			}
			else if ((x==350 || x==451) && (y>290 && y<295)) {
				pixmem32 = (Uint32*) (*surf)->pixels + y*(*surf)->w + x;
				*pixmem32 = SDL_MapRGBA((*surf)->format, 255, 255, 255, 0);
			}
			else {
				pixmem32 = (Uint32*) (*surf)->pixels + y*(*surf)->w + x;
				*pixmem32 = SDL_MapRGBA((*surf)->format, 0, 0, 0, 200);
			}
		}
	}
	for (y=291;y<295;y++) {
		for (x=350;x<=350+(int)((double)i/(double)imax*(451.0-350.0));x++) {
			pixmem32 = (Uint32*) (*surf)->pixels + y*(*surf)->w + x;
			*pixmem32 = SDL_MapRGBA((*surf)->format, 255, 255, 255, 0);
		}
	}
	surf_tex = SDL_CreateTextureFromSurface(renderer, (*surf));
	loading_tex = renderText(text, FontFile, white, 14, renderer);

	SDL_RenderClear(renderer);
	ApplySurface(0, 0, bgtex, renderer, NULL);
	ApplySurface(0, 0, surf_tex, renderer, NULL);
	renderTexture(loading_tex, renderer, 275, 260);
	SDL_RenderPresent(renderer);

	SDL_DestroyTexture(loading_tex);
	SDL_DestroyTexture(surf_tex);

	return 0;
}

#endif /* PLOT_H_ */
