#pragma once

#include <X11/X.h>
#include <X11/Xlib.h>
#include <bit>
#include <concepts>
#include <cstring>
#include <vector>

#include "check.hpp"
#include "types.hpp"

struct Rect
{
    i32 left, top;
    u32 width, height;
};

class XWindow
{
  public:
    static XWindow create(size_t width, size_t height)
    {
        XWindow window;

        window.display = XOpenDisplay(NULL);
        if (!window.display)
        {
            std::cerr << "Unable to open Display" << std::endl;
            std::exit(1);
        }

        window.screen = DefaultScreen(window.display);

        window.window = XCreateSimpleWindow(window.display, RootWindow(window.display, window.screen), 0, 0, (u32)width,
                                            (u32)height, 0, BlackPixel(window.display, window.screen),
                                            WhitePixel(window.display, window.screen));

        XStoreName(window.display, window.window, "Window");
        XSelectInput(window.display, window.window, ExposureMask | KeyPressMask);

        window.gc = XCreateGC(window.display, window.window, 0, NULL);
        if (!window.gc)
        {
            std::cerr << "Unable to create GC" << std::endl;
            std::exit(1);
        }

        XSetForeground(window.display, window.gc, BlackPixel(window.display, window.screen));

        XMapWindow(window.display, window.window);

        window.backing_buffer.resize(width * height, 0u);
        window.image =
            XCreateImage(window.display, DefaultVisual(window.display, window.screen),
                         (u32)DefaultDepth(window.display, window.screen), ZPixmap, 0,
                         reinterpret_cast<char *>(window.backing_buffer.data()), (u32)width, (u32)height, 32, 0);

        window.pixmap = XCreatePixmap(window.display, window.window, (u32)width, (u32)height,
                                      (u32)DefaultDepth(window.display, window.screen));

        return window;
    }

    ~XWindow()
    {
        XFreePixmap(display, pixmap);
        XFreeGC(display, gc);
        XDestroyWindow(display, window);
        XCloseDisplay(display);
    }

    std::vector<XEvent> poll_events()
    {
        XEvent event;
        std::vector<XEvent> events;
        size_t event_count = (size_t)XPending(display);
        events.reserve(event_count);
        for (size_t i = 0; i < event_count; ++i)
        {
            XNextEvent(display, &event);
            events.push_back(event);
        }
        return events;
    }

    Rect active_area() const
    {
        XWindowAttributes attrs;
        XGetWindowAttributes(display, window, &attrs);
        i32 border = attrs.border_width;

        CHECK(attrs.width >= border);
        CHECK(attrs.height >= border);

        u32 width = (u32)(attrs.width - border);
        u32 height = (u32)(attrs.height - border);
        return {border, border, width, height};
    }

    friend class DrawFrame;

  public:
    bool should_close = false;

  private:
    Display *display;
    int screen;
    Window window;
    GC gc;
    XImage *image;
    Pixmap pixmap;
    std::vector<u32> backing_buffer{};
};

class DrawFrame
{
  public:
    DrawFrame(XWindow &window_) : window(window_)
    {
        auto display = window.display;
        auto gc = window.gc;
        auto screen = window.screen;
        auto rect = window.active_area();
        XSetForeground(display, gc, WhitePixel(display, screen));
        XFillRectangle(display, window.pixmap, gc, rect.left, rect.top, rect.width, rect.height);
        XSetForeground(display, gc, BlackPixel(display, screen));
    }

    ~DrawFrame()
    {
        auto display = window.display;
        auto gc = window.gc;
        auto rect = window.active_area();
        auto pixmap = window.pixmap;
        auto image = window.image;

        XFlush(display);
        XPutImage(display, pixmap, gc, image, 0, 0, 0, 0, (u32)image->width, (u32)image->height);
        XCopyArea(display, pixmap, window.window, gc, 0, 0, rect.width, rect.height, 0, 0);
    }

    template <std::floating_point T> void blit(const Image<RGB<T>> &image)
    {
        using namespace std;

        CHECK(image.width() == (size_t)window.image->width);
        CHECK(image.height() == (size_t)window.image->height);
        const auto rmask = (u32)window.image->red_mask;
        const auto gmask = (u32)window.image->green_mask;
        const auto bmask = (u32)window.image->blue_mask;

        const auto rshift = countr_zero(rmask);
        const auto gshift = countr_zero(gmask);
        const auto bshift = countr_zero(bmask);

        const auto pack_pixel = [&](const RGB<T> &color) {
            const auto r = std::clamp(color.r * T{255}, T{0}, T{255});
            const auto g = std::clamp(color.g * T{255}, T{0}, T{255});
            const auto b = std::clamp(color.b * T{255}, T{0}, T{255});
            return ((u32)r << rshift & rmask) | ((u32)g << gshift & gmask) | ((u32)b << bshift & bmask);
        };

        u32 *dst = reinterpret_cast<u32 *>(window.image->data);
        transform(begin(image), end(image), dst, pack_pixel);
    }

  private:
    XWindow &window;
};
