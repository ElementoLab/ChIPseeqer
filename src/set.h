#ifndef SET_H
#define SET_H

#define ALLOC_BLOCK	32

/* opaque type for the set */
typedef struct set {
	char **arr;
	size_t size;
	size_t capacity;
} set_t;

#endif /* SET_H */
