#ifndef HEAP_H_
#define HEAP_H_

/* Last updated: Dec 31, 2015
 * Copyright (c) 2015 Antti Hakkinen
 */

#ifdef __cplusplus
extern "C" {
#endif

#include <assert.h>

#define HEAP_DEFINE(_PREFIX, _TYPE, _SIZE, _LESS, _COOKIE_TYPE) \
	\
	static \
	_SIZE \
	_PREFIX ## _heapdn(_TYPE *x, _SIZE n, _SIZE i, _TYPE v, _COOKIE_TYPE cookie) \
	{ \
		_SIZE j; \
		j = 2*i; \
		while (!(j+1 > n)) \
		{ \
			if (_LESS(x[j-1], x[j+1-1], cookie)) \
				++j; \
			if (!_LESS(v, x[j-1], cookie)) \
				goto done; \
			x[i-1] = x[j-1]; \
			i = j; \
			j = 2*i; \
		} \
		if (!(j > n)) \
		{ \
			if (!_LESS(v, x[j-1], cookie)) \
				goto done; \
			x[i-1] = x[j-1]; \
			i = j; \
		} \
	done: \
		x[i-1] = v; \
		return i; \
	} \
	\
	static \
	_SIZE \
	_PREFIX ## _heapup(_TYPE *x, _SIZE j, _TYPE v, _COOKIE_TYPE cookie) \
	{ \
		_SIZE i; \
		while (j > 1) \
		{ \
			i = j/2; \
			if (!_LESS(x[i-1], v, cookie)) \
				break; \
			x[j-1] = x[i-1]; \
			j = i; \
		} \
		x[j-1] = v; \
		return j; \
	} \
	\
	static \
	_SIZE \
	_PREFIX ## _heapholedn(_TYPE *x, _SIZE n, _SIZE i, _COOKIE_TYPE cookie) \
	{ \
		_SIZE j; \
		j = 2*i; \
		while (!(j+1 > n)) \
		{ \
			if (_LESS(x[j-1], x[j+1-1], cookie)) \
				++j; \
			x[i-1] = x[j-1]; \
			i = j; \
			j = 2*i; \
		} \
		if (!(j > n)) \
		{ \
			x[i-1] = x[j-1]; \
			i = j; \
		} \
		return i; \
	} \
	\
	static \
	void \
	_PREFIX ## _make_heap(_TYPE *x, _SIZE n, _COOKIE_TYPE cookie) \
	{ \
		_SIZE t; \
		assert(n > 0); \
		for (t = n/2; t > 0; --t) \
			_PREFIX ## _heapdn(x, n, t, x[t-1], cookie); \
	} \
	\
	static \
	void \
	_PREFIX ## _push_heap(_TYPE *x, _SIZE n, _COOKIE_TYPE cookie) \
	{ \
		assert(n > 0); \
		_PREFIX ## _heapup(x, n, x[n-1], cookie); \
	} \
	\
	static \
	void \
	_PREFIX ## _pop_heap(_TYPE *x, _SIZE n, _COOKIE_TYPE cookie) \
	{ \
		_TYPE v; \
		_SIZE j; \
		assert(n > 0); \
		v = x[1-1]; \
		j = _PREFIX ## _heapholedn(x, n, 1, cookie); \
		if (j < n) \
			_PREFIX ## _heapup(x, j, x[n-1], cookie); \
		x[n-1] = v; \
	} \
	\
	static \
	void \
	_PREFIX ## _sort_heap(_TYPE *x, _SIZE n, _COOKIE_TYPE cookie) \
	{ \
		_SIZE t; \
		assert(n > 0); \
		for (t = n; t > 1; --t) \
			_PREFIX ## _pop_heap(x, t, cookie); \
	} \

#ifdef __cplusplus
}
#endif

#endif /* HEAP_H_ */
