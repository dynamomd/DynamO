/*    dynamo:- Event driven molecular dynamics simulator 
 *    http://www.dynamomd.org
 *    Copyright (C) 2009  Marcus N Campbell Bannerman <m.bannerman@gmail.com>
 *
 *    This program is free software: you can redistribute it and/or
 *    modify it under the terms of the GNU General Public License
 *    version 3 as published by the Free Software Foundation.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

//Here we have the correct order of GL includes
#include <magnet/GL/context.hpp>
#include <SDL2/SDL.h>

namespace magnet
{
    namespace GL
    {
        namespace SDL
        {
            namespace detail
            {
                void checkSDLError()
                {
#ifndef NDEBUG
                    const char *error = SDL_GetError();
                    if (*error != '\0')
                        M_throw() << "SDL Error " << error;
#endif
                }
            }

            class Window
            {
            protected:
                SDL_Window* _handle;
                SDL_GLContext _context;
                Sint32 _w;
                Sint32 _h;

            public:
                Window() {}

                void init(std::string name, int w = 640, int h = 480, Uint32 flags = SDL_WINDOW_RESIZABLE, SDL_GLContext context = nullptr)
                {
                    _context = context;
                    _w = w;
                    _h = h;
                    _handle = SDL_CreateWindow(name.c_str(), SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED,
                                               _w, _h, flags | SDL_WINDOW_OPENGL);

                    if (!_handle) /* Die if creation failed */
                    {
                        detail::checkSDLError();
                        M_throw() << "Unable to create window";
                    }
                }

                SDL_Window* getSDLHandle() const { return _handle; }

                void setContext(SDL_GLContext context)
                {
                    _context = context;
                }

                SDL_GLContext getContext() const { return _context; }

                void makeCurrent()
                {
                    if (SDL_GL_MakeCurrent(_handle, _context))
                        detail::checkSDLError();
                }

                void deinit()
                {
                    SDL_DestroyWindow(_handle);
                }
            };

            class Engine
            {
            protected:
                volatile bool _running = false;
                /* Our opengl context handle */
                Window _mainWindow;
                magnet::GL::Context::ContextPtr _context;
                std::thread _renderer;

            public:
                Engine() {}

                void init(std::string windowName = "")
                {
                    if (SDL_Init(SDL_INIT_VIDEO) < 0) /* Initialize SDL's Video subsystem */
                    {
                        SDL_Quit();
                        M_throw() << "Unable to initialize SDL"; /* Or die on error */
                    }
                    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION, 3);
                    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION, 2);

                    /* Create our window centered at 512x512 resolution */
                    _mainWindow.init("Main window", 640, 480, SDL_WINDOW_RESIZABLE);
                    
                    {
                    /* Create our opengl context and attach it to our window */
                    SDL_GLContext glsdlcontext =  SDL_GL_CreateContext(_mainWindow.getSDLHandle());

                    if (!glsdlcontext)
                        detail::checkSDLError();
                    _mainWindow.setContext(glsdlcontext);
                    }

                    _running = true;
                    _context = magnet::GL::Context::getContext();
                }

                void deinit()
                {
                    /* Delete our opengl context, destroy our window, and shutdown SDL */
                    SDL_GL_DeleteContext(_mainWindow.getContext());
                    _mainWindow.deinit();
                    SDL_Quit();
                }

                void tick()
                {
                    SDL_Event event;
                    while (SDL_PollEvent(&event))
                        switch (event.type)
                        {
                        case SDL_QUIT:
                            _running = false;
                            break;
                        default:
                            M_throw() << "Unhandled SDL event " << event.type;
                        }

                    _context->tick();
                }

                void go_render()
                {
                    init();
                    while (_running)
                        tick();
                    deinit();
                }

                void launch_render_thread()
                {
                    _renderer = std::thread(&Engine::go_render, this);
                    return;
                }

                magnet::GL::Context::ContextPtr getContext()
                {
                    return _context;
                }

                std::thread &getRenderthread()
                {
                    return _renderer;
                }
            };
        }
    }
}
