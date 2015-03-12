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
int scanNumber(char (*nb)[10], int i);
int File2tex(char* title, SDL_Texture** tex, char* path);
int File2surf(char* title, SDL_Surface** surf, char* path);

SDL_Window* window;
SDL_Renderer* renderer;

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

// Get a pixel

Uint32 get_pixel32( SDL_Surface *surface, int x, int y ) {
    //Convert the pixels to 32 bit
    Uint32 *pixels = (Uint32 *)surface->pixels;

    //Get the requested pixel
    return pixels[ ( y * surface->w ) + x ];
}

// Put a pixel

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

// Convert number to string, right-aligned
int scanNumber(char (*nb)[10], int i) {
	if (i<100) sprintf((*nb), "   %d", i);
	else if (i<1000) sprintf((*nb), "  %d", i);
	else sprintf((*nb), "%d", i);
	return 0;
}

// Open image file into a texture
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

// Open image file into a surface
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

#endif /* PLOT_H_ */
