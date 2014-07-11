#include "include/generaldefinitions.h"

bool isInteger(double value)
{
    return (value == floor(value)); // Fast, to int precision
    // return (abs(value - round(value)) < 1e-9); // Slow, specified precision
}

// Returns a random integer in the range [min, max]
int randomInteger(int min, int max)
{
    static bool randomSeedSet = false;
    if (!randomSeedSet)
    {
        srand(time(NULL)); // Random seed
        randomSeedSet = true;
    }
    return rand() % (max - min + 1) + min;
}
