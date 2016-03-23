/*                   _           
 *   _ __ ___   ___ | | ___  ___ 
 *  | '_ ` _ \ / _ \| |/ _ \/ __|
 *  | | | | | | (_) | |  __/ (__ 
 *  |_| |_| |_|\___/|_|\___|\___| - Molecular Dynamics Framework
 *          
 *  Copyright (C) 2016  Carlo Del Don  (deldonc@student.ethz.ch)
 *                      Michel Breyer  (mbreyer@student.ethz.ch)
 *                      Florian Frei   (flofrei@student.ethz.ch)
 *                      Fabian Thuring (thfabian@student.ethz.ch)
 * 
 *  This file is distributed under the MIT Open Source License. 
 *  See LICENSE.txt for details.
 */

#ifndef MOLEC_CONFIG_H
#define MOLEC_CONFIG_H

#if defined(__clang__)
#define MOLEC_COMPILER_CLANG 1
#endif

#if defined(__ICC) || defined(__INTEL_COMPILER)
#define MOLEC_COMPILER_INTEL 1
#endif 

#if defined(__GNUC__) || defined(__GNUG__)
#define MOLEC_COMPILER_GNU 1
#endif

#if defined(__PGI)
#define MOLEC_COMPILER_PGI 1
#endif

#if defined(_MSC_VER)
#define MOLEC_COMPILER_MSVC 1
#endif

#if defined(_MSC_VER) || defined(_WIN32) || defined(_WIN64)
#define MOLEC_PLATFORM_WINDOWS 1
#elif defined(__linux__ ) || defined(__linux)
#define MOLEC_PLATFORM_LINUX 1
#elif defined(__APPLE__)
#define MOLEC_PLATFORM_APPLE 1
#endif

#if defined (__unix__) || defined(MOLEC_PLATFORM_APPLE)
#define MOLEC_PLATFORM_POSIX 1
#endif

// ASM
#if defined(MOLEC_COMPILER_GNU)
#define MOLEC_ASM __asm__
#else
#define MOLEC_ASM __asm__ 
#endif

// VOLATILE
#if defined(MOLEC_COMPILER_GNU)
#define MOLEC_VOLATILE __volatile__
#else
#define MOLEC_VOLATILE volatile
#endif

// NORETURN
#if defined(MOLEC_COMPILER_GNU)
#define MOLEC_NORETURN __attribute__((noreturn))
#elif defined(MOLEC_COMPILER_MSVC)
#define MOLEC_NORETURN __declspec(noreturn)
#else
#define MOLEC_NORETURN
#endif

// INLINE
#if defined(MOLEC_COMPILER_GNU)
#define MOLEC_INLINE inline __attribute__((always_inline))
#elif defined(MOLEC_COMPILER_MSVC)
#define MOLEC_INLINE __forceinline
#else
#define MOLEC_INLINE inline
#endif

// NOINLINE
#if defined(MOLEC_COMPILER_GNU)
#define MOLEC_NOINLINE __attribute__((noinline))
#elif defined(MOLEC_COMPILER_MSVC)
#define MOLEC_NOINLINE __declspec(noinline)
#else
#define MOLEC_NOINLINE
#endif

// Disable some warnings on Windows
#ifdef MOLEC_PLATFORM_WINDOWS
#define _CRT_SECURE_NO_WARNINGS

#pragma warning(disable : 4267) // conversion from 'size_t' to 'int'
#pragma warning(disable : 4244) // conversion from 'double' to 'int'
#pragma warning(disable : 4305) // truncation from 'double' to 'float'
#endif

#endif
