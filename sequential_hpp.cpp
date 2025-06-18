#include <iostream>
#include <cstdint> // für uint8_t
#include <bitset>  // für std::bitset

void printBits(uint8_t value)
{
    std::bitset<8> bits(value);
    std::cout << bits;
}

uint8_t collision(uint8_t current_cell)
{
    if (current_cell == 0b00000101)
    {
        current_cell = 0b00001010; // Setze a auf 0b00001010
    }
    else if (current_cell == 0b00001010)
    {
        current_cell = 0b00000101; // Setze a auf 0b00000101
    }
    return current_cell; // Rückgabe des aktualisierten Wertes
}

int main()
{                           // ---------0bxxxknosw------------------- k=wall bitt
    uint8_t a = 0b00000101; // Binärliterale ab C++14 erlaubt
    uint8_t b = 0b00001110; // Binärliterale ab C++14 erlaubt
    std::cout << "Start-Wert a: ";
    printBits(a);                                         // Ausgabe der Bits von a
    std::cout << " (dezimal: " << +a << ")" << std::endl; // + damit es als Zahl ausgegeben wird
    std::cout << "Start-Wert b: ";
    printBits(b);                                         // Ausgabe der Bits von b
    std::cout << " (dezimal: " << +b << ")" << std::endl; // + damit es als Zahl ausgegeben wird

    a = collision(a); // Aufruf der Funktion collision mit a
    b = collision(b); // Aufruf der Funktion collision mit b

    std::cout << "End-Wert a: ";
    printBits(a);                                         // Ausgabe der Bits von a
    std::cout << " (dezimal: " << +a << ")" << std::endl; // + damit es als Zahl ausgegeben wird
    std::cout << "End-Wert b: ";
    printBits(b);                                         // Ausgabe der Bits von b
    std::cout << " (dezimal: " << +b << ")" << std::endl; // + damit es als Zahl ausgegeben wird
    return 0;
}