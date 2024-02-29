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
    0x107c,
    0x9237,
    0xdc9d,
    0x6422,
    0x7c21,
    0x2033,
    0x5572,
    0x85be,
    0xc5af,
    0x3f16,
    0xa257,
    0x7bca,
    0x34f9,
    0x7334,
    0x143f,
    0x51de,
    0xfd3e,
    0x2ee6,
    0x3ef8
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
    
    static uint32_t seedValue = 0XAABBCCDD;
    uint32_t key = CalKey(seedValue);
    // cout << key << endl;
    // printf("%X", key);
    // printf("%X", 0x8c673545);
   

    return 0;
}
