/******************************************************************************

                              Online C++ Compiler.
               Code, Compile, Run and Debug C++ program online.
Write your code in this editor and press "Run" button to compile and execute it.

*******************************************************************************/

#include <iostream>

using namespace std;



#define keyLen 0x13u

#define ODDEVENCHECK15 0x8000u
#define ODDEVENCHECK4 0x0010u
#define ODDEVENCHECK2 0x0004u
#define ODDEVENCHECK1 0x0002u

#define BIT15 15u
#define BIT4 4u
#define BIT2 2u
#define BIT1 1u

static const uint16_t keydata[0x13] = 
{
    0xb5d5,
    0x1741,
    0xfe95,
    0x4ad5,
    0x9497,
    0xddc5,
    0x3761,
    0x4b0d,
    0x5f5d,
    0x6b23,
    0xcfa9,
    0xc8e7,
    0xf209,
    0x4fff,
    0xff2b,
    0x4801,
    0x8891,
    0x48b9,
    0x3607
};


static uint32_t CalKey(uint32_t seedValue);
static inline uint16_t cc16_from_u8(uint16_t high, uint16_t low);
static uint16_t Fun_XOR(uint16_t Data1, uint16_t Data2);
static uint16_t Fun_SwapData(uint16_t Data);
static inline uint16_t cc_8high(uint32_t value);
static inline uint16_t cc_8low(uint32_t value);
static inline uint16_t cc_16high(uint32_t value);
static inline uint16_t cc_16low(uint32_t value);
static inline uint16_t cc32_from_u16(uint16_t high, uint16_t low);


static inline uint16_t cc16_from_u8(uint16_t high, uint16_t low)
{
    return (((uint16_t) high) * (255 + 1)) + low;
}

static uint16_t Fun_XOR(uint16_t Data1, uint16_t Data2)
{
    return (Data1 ^ Data2);
}

static uint16_t Fun_SwapData(uint16_t Data)
{
    uint8_t HighByte;
    uint8_t LowerByte;
    uint16_t return_U;
    HighByte = cc_8high(Data);
    LowerByte = cc_8low(Data);
    return_U = cc16_from_u8(LowerByte, HighByte);
    return return_U;
}

static inline uint16_t cc_16high(uint32_t value)
{
    return value / (((uint32_t) 65535 + 1));
}

static inline uint16_t cc_16low(uint32_t value)
{
    return value % (((uint32_t) 65535 + 1));
}

static inline uint16_t cc_8high(uint32_t value)
{
    return value / (255 + 1);
}

static inline uint16_t cc_8low(uint32_t value)
{
    return value % (255 + 1);
}

static inline uint16_t cc32_from_u16(uint16_t high, uint16_t low)
{
    return (((uint32_t) high) * ((uint32_t) 65536)) + low;
}

static uint32_t CalKey(uint32_t seedValue)
{
    uint8_t ilterCount;
    uint8_t oddEvenCheck;
    uint16_t l_InputHighOrderSeed;
    uint16_t l_InputLowOrderSeed;
    uint16_t l_OutputLowOrderSeed;
    
    l_InputHighOrderSeed = cc_16high(seedValue);
    l_InputLowOrderSeed = cc_16low(seedValue);

    // printf("%X", l_InputHighOrderSeed);
    // printf("%X", l_InputLowOrderSeed);
    
    for (ilterCount = 0u; ilterCount < keyLen;)
    {
        
        l_OutputLowOrderSeed = Fun_XOR(l_InputLowOrderSeed, keydata[ilterCount]);
        oddEvenCheck = ((l_OutputLowOrderSeed & ODDEVENCHECK15) >> BIT15) + ((l_OutputLowOrderSeed & ODDEVENCHECK4) >> BIT4) + ((l_OutputLowOrderSeed & ODDEVENCHECK2) >> BIT2) + ((l_OutputLowOrderSeed & ODDEVENCHECK1) >> BIT1);
        l_OutputLowOrderSeed = l_OutputLowOrderSeed * 2u;
        if (oddEvenCheck % 2)
        {
            l_OutputLowOrderSeed = l_OutputLowOrderSeed + 1u;
        } 
        else 
        {
            
        }
        
        // printf("%X", l_OutputLowOrderSeed);
        // printf("\n");

        l_OutputLowOrderSeed = Fun_SwapData(l_OutputLowOrderSeed);
        l_OutputLowOrderSeed = Fun_XOR(l_OutputLowOrderSeed, keydata[ilterCount]);
        
        l_InputLowOrderSeed = l_InputHighOrderSeed;
        l_InputHighOrderSeed = Fun_XOR(l_InputLowOrderSeed, l_OutputLowOrderSeed);
        ilterCount = ilterCount + 1;
    }

    printf("%X", l_InputHighOrderSeed);
    printf("\n");
    printf("%X", l_InputLowOrderSeed);
    printf("\n");
    
    
    return cc32_from_u16(l_InputHighOrderSeed, l_InputLowOrderSeed);
}


int main()
{
    // cout<<"Hello World";
    
    static uint32_t seedValue = 0X11221122;
    uint32_t key = CalKey(seedValue);
    // cout << key << endl;
    // printf("%X", key);
    // printf("%X", 0x8c673545);
   

    return 0;
}
