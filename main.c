#include <X11/Xlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <unistd.h>
#include <poll.h>
#include <time.h>

// #define MAX_FPS 5

#define WIDTH  800
#define HEIGHT 800

#define WINDOW_BORDER_WIDTH 1


// ** Pixels **
typedef uint32_t Pixel;

#define PIXEL_DEPTH 24

// hard coded for Little-Endian
#define PIXEL_RGB(r, g, b) (Pixel) ( \
	(r << 8 * 2) |					 \
	(g << 8    ) |					 \
	(b		   )  					 \
)
#define PIXEL_R(p) ((p >> 8 * 2) & 0xff)
#define PIXEL_G(p) ((p >> 8    ) & 0xff)
#define PIXEL_B(p) ((p		   ) & 0xff)

// ** Pixel Buffer **
Pixel pixel_buffer[WIDTH * HEIGHT] = {0};
	


// ** Vectors **
typedef struct V3 {
	float x;
	float y;
	float z;
} V3;

float dot_product(V3 a, V3 b) {
	return (
		a.x * b.x +
	   	a.y * b.y + 
		a.z * b.z
	);
}

float norm(V3 a) {
	return sqrt(
		a.x * a.x +
		a.y * a.y +
		a.z * a.z
	);
}

V3 normalize(V3 a) {
	float n = norm(a);

	return (V3) {
		a.x / n,
		a.y / n,
		a.z / n
	};
}

V3 scalar_product(V3 a, float b) {
	return (V3) {
		b * a.x,
		b * a.y,
		b * a.z
	};
}

V3 vector_sum(V3 a, V3 b) {
	return (V3) {
		a.x + b.x,
		a.y + b.y,
		a.z + b.z
	};
}
V3 vector_sub(V3 a, V3 b) {
	return (V3) {
		a.x - b.x,
		a.y - b.y,
		a.z - b.z
	};
}


// ** World **
typedef struct Sphere {
	V3 position; 	// center point

	float radius;
} Sphere;

typedef struct LightSource {
	V3 position;

	Pixel colour;
	float intensity;
} LightSource;

typedef struct Camera {
	V3 position; 	// eye (ray origin) position
							   
	// the plane has its center on the local x-axis 
	// and it's orhogonal to it
	float near_plane_distance; // relative to eye
	float near_plane_width;
	float near_plane_height;
} Camera;

int intersections_ray_sphere(Sphere sphere, V3 ray_origin, V3 ray_direction, V3* intersections) {
	assert(intersections);

	// initialize results
	intersections[0] = (V3){
		-INFINITY,
		-INFINITY,
		-INFINITY
	};
	intersections[1] = (V3){
		-INFINITY,
		-INFINITY,
		-INFINITY
	};
	int ret = 0;

	// normalize ray direction vector
	ray_direction = normalize(ray_direction);

	V3 _oc = vector_sub(ray_origin, sphere.position);

	float _alpha = dot_product(ray_direction, _oc);
	float _beta  = (dot_product(_oc, _oc) - sphere.radius * sphere.radius);
	float _delta = _alpha * _alpha - _beta;

	if(_delta < 0) {
		return 0; // there are no intersections
	}

	float t1 = -_alpha + sqrt(_delta);
	float t2 = -_alpha - sqrt(_delta);

	// calculate intersection points
	if(t1 >= 0) {
		intersections[0] = vector_sum(ray_origin, scalar_product(ray_direction, t1));
		ret++;
	}
	if(t2 >= 0) {
		intersections[1] = vector_sum(ray_origin, scalar_product(ray_direction, t2));
		ret++;
	}

	return ret;
}

// ** Stupid Simple render **
int render_sphere(Camera camera, Sphere sphere, LightSource light, Pixel* out_buffer, size_t width, size_t height) {
	assert(out_buffer);
	assert(width  > 0);
	assert(height > 0);

	// TODO: world coordinates
	// for now we just assert that the camera position is (0, 0, 0), because we haven't figured out actual world coordinates yet.
	assert(!(
		camera.position.x ||
		camera.position.y ||
		camera.position.z
	));


	// define center of near plane
	V3 near_plane_center = {
		camera.near_plane_distance,
		0,
		0
	};

	// for every pixel (i, j)
	for(size_t i = 0; i < height; i++) {
		for(size_t j = 0; j < width; j++) {
			out_buffer[i * width + j] = PIXEL_RGB(0x00, 0x00, 0x00);

			// figure out where (i, j) is on the near plane
			float near_plane_i = (float)i / height * camera.near_plane_height;
			float near_plane_j = (float)j / width  * camera.near_plane_width;
			float near_plane_x = near_plane_j - camera.near_plane_width  / 2.0;
			float near_plane_y = near_plane_i - camera.near_plane_height / 2.0;

			// ray origin
			V3 ray_origin = camera.position;

			// calculate direction vector of the ray
			V3 ray_direction = {
				near_plane_center.x,
				near_plane_x,
				near_plane_y
			};

			// calculate intersections
			V3 intersections[2];
			
			if(!intersections_ray_sphere(sphere, ray_origin, ray_direction, intersections)) {
				continue;
			}

			// we only care about the closer intersection. Only think about the x value for now.
			V3 closer_intersection = intersections[0];
			if(intersections[0].x > intersections[1].x) {
				closer_intersection = intersections[1];
			}
			assert(closer_intersection.x != -INFINITY);

			// calculate normal of the sphere at the intersection
			V3 sphere_intersection_normal = normalize(vector_sub(closer_intersection, sphere.position));
			
			// find intersection -> light source vector
			ray_direction = vector_sub(light.position, closer_intersection);
			ray_origin = closer_intersection;

			// check for occlusions (shadows)
			// intersections_ray_sphere(sphere, ray_origin, ray_direction, intersections);

			// // we check wether there is an intersection that is distant enough from the origin
			// if(intersections[0].x != -INFINITY && norm(vector_sub(intersections[0], ray_origin)) >= 0.01) {
			// 	continue;
			// }
			// if(intersections[1].x != -INFINITY && norm(vector_sub(intersections[1], ray_origin)) >= 0.01) {
			// 	continue;
			// }

			
			// now we can do the dot product of the surface normal with the direction of the light
			// making sure to normalize the direction of the light
			float reflectance = dot_product(sphere_intersection_normal, normalize(ray_direction));

			// assert(reflectance <= 1.0);
			if(reflectance < 0.0) {
				// the surface is facing away from the light source
				continue;
			}

			// assert(light.intensity <= 1.0 && light.intensity >= 0.0);
			Pixel p = PIXEL_RGB(
				(Pixel) (PIXEL_R(light.colour) * reflectance * light.intensity),
				(Pixel) (PIXEL_G(light.colour) * reflectance * light.intensity),
				(Pixel) (PIXEL_B(light.colour) * reflectance * light.intensity)
			);

			out_buffer[i * width + j] = p;
		}
	}
	

	return 0;
}

int main() {
	// ** X display and window **
    Display *display;
    Window window;
	int screen;
	
    display = XOpenDisplay(NULL);
	if(!display) {
		printf("error: could not open X11 display\n");
		return 1;
	}

    screen = DefaultScreen(display);

    window = XCreateSimpleWindow(
		display, RootWindow(display, screen),
	   	0, 0, 
		WIDTH, HEIGHT, 
		WINDOW_BORDER_WIDTH,
		BlackPixel(display, screen), WhitePixel(display, screen)
	);

	// ** X events **
    XSelectInput(display, window, ExposureMask | KeyPressMask);
	

	// ** X graphics context **
	GC gc = XCreateGC(display, window, 0, NULL);
	

	// ** X image **
	XImage* image = XCreateImage(
		display, 
		DefaultVisual(display, screen),
	   	PIXEL_DEPTH, 
		ZPixmap, 
		0,
        (char*) pixel_buffer, 
		WIDTH, HEIGHT, 
		8 * sizeof(Pixel), 
		0
	);
	if(!image) {
		printf("error: could not create XImage.\n");
		return 1;
	}


	// create some objects
	Camera cute_camera = {
		.position = (V3){0, 0, 0},
		.near_plane_distance = 2,
		.near_plane_width  = 1.8, // (float)WIDTH  / 1000.0,
		.near_plane_height = 1.8  // (float)HEIGHT / 1000.0
	};

	// let's think locally to the camera for now. Every position is relative to the camera eye.
	Sphere sphere = {
		.position = {20.0, 0.0, 0.0},
		.radius   = 3.0
	};
	LightSource light = {
		.position  = {0.0, 0.0, 20.0},
		.colour = PIXEL_RGB(0x65, 0x32, 0x26),
		.intensity = 1.0
	};


	// ** Event loop **
    XMapWindow(display, window);

	// register interest in the delete window message
	Atom wm_delete_window = XInternAtom(display, "WM_DELETE_WINDOW", False);
	XSetWMProtocols(display, window, &wm_delete_window, 1);

	// find x11 events file descriptor
	struct pollfd xfd = {
		.fd = ConnectionNumber(display),
		.events = POLLIN,
	};

	// time setup
	struct timespec lt, ct;
	double dt, game_time;

	clock_gettime(CLOCK_MONOTONIC, &lt);

	int should_quit = 0;
	while (!should_quit) {
		// elapsed time
		clock_gettime(CLOCK_MONOTONIC, &ct);
		dt = (ct.tv_sec - lt.tv_sec) + (ct.tv_nsec - lt.tv_nsec) / 1e9;
		lt = ct;

		game_time += dt;
	
		printf("FPS : %f\n", 1.0 / dt);
			
		while(XPending(display)) {
			XEvent event;
			XNextEvent(display, &event);

			switch (event.type) {
				case ClientMessage:
					if ((Atom)event.xclient.data.l[0] == wm_delete_window)
						should_quit = 1;
					break;
			}
		}

		int num_ready = poll(&xfd, 1, 15);
		if(num_ready < 0) {
			perror("poll()");
			break;
		}

		light.position.y = 20.0 * sin(game_time * 2.0 * M_PI * 0.2);
		light.position.z = 20.0 * cos(game_time * 2.0 * M_PI * 0.2);
		// light.colour = PIXEL_RGB(
		// 	(Pixel) ((1.0 + sin(game_time * 2.0 * M_PI * 0.2)) * 0x7f),
		// 	(Pixel) ((1.0 + cos(game_time * 2.0 * M_PI * 0.2)) * 0x7f),
		// 	(Pixel) ((1.0 + sin(game_time * 2.0 * M_PI * 0.2)) * 0x7f)
		// );
		render_sphere(cute_camera, sphere, light, pixel_buffer, WIDTH, HEIGHT);

		XPutImage(display, window, gc, image, 0, 0, 0, 0, WIDTH, HEIGHT);
		XFlush(display);

		#ifdef MAX_FPS
		{
			double time_to_sleep = 1.0 / MAX_FPS - dt;
			struct timespec ts;
			if(time_to_sleep > 0) {
				ts.tv_sec  = (time_t)time_to_sleep;
				ts.tv_nsec = (long)((time_to_sleep - ts.tv_sec) * 1e9);
				nanosleep(&ts, NULL);
			}
		}
		#endif
	}


    XDestroyWindow(display, window);
    XCloseDisplay(display);
    return 0;
}
