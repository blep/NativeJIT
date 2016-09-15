#include <chrono>
#include <iostream>
#include <fstream>
#include <math.h>       // For value of e.
#include <stdexcept>
#include <string>
#include <sstream>

#include "NativeJIT/CodeGen/ExecutionBuffer.h"
#include "NativeJIT/CodeGen/FunctionBuffer.h"
#include "NativeJIT/Function.h"
#include "Temporary/Allocator.h"


using NativeJIT::Allocator;
using NativeJIT::ExecutionBuffer;
using NativeJIT::Function;
using NativeJIT::FunctionBuffer;


namespace Examples
{

    struct PerfStat
    {
        long long parseNano_;
        long long compileNano_;
        long long evalNano_;

        void print() 
        {
            std::cout << "Parse time: " << parseNano_ << " ns, Compile time: " << compileNano_ << " ns, Eval time: " << evalNano_ << " ns" << std::endl;
        }
    };

#if 1 /* idea to avoid templated parsing function */

    template<NativeJIT::JccType JCC, typename T, typename U>
    NativeJIT::Node<T> &makeIf( Function<T> &expr, NativeJIT::Node<U> &cmpLhs, NativeJIT::Node<U> &cmpRhs,
                                NativeJIT::Node<T> &trueValue, NativeJIT::Node<T> &falseValue )
    {
        auto &cmp = expr.Compare<JCC>( cmpLhs, cmpRhs );
        return expr.If( cmp, trueValue, falseValue );
    }

    template<typename T, typename U>
    NativeJIT::Node<T> &If( Function<T> &expr, NativeJIT::Node<U> &cmpLhs, NativeJIT::JccType cmpOp, NativeJIT::Node<U> &cmpRhs,
                            NativeJIT::Node<T> &trueValue, NativeJIT::Node<T> &falseValue )
    {
        switch ( cmpOp )
        {
        case NativeJIT::JccType::JO:
            return makeIf<NativeJIT::JccType::JO>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JNO:
            return makeIf<NativeJIT::JccType::JNO>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JB:
            return makeIf<NativeJIT::JccType::JB>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JAE:
            return makeIf<NativeJIT::JccType::JAE>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JE:
            return makeIf<NativeJIT::JccType::JE>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JNE:
            return makeIf<NativeJIT::JccType::JNE>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JBE:
            return makeIf<NativeJIT::JccType::JBE>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JA:
            return makeIf<NativeJIT::JccType::JA>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JS:
            return makeIf<NativeJIT::JccType::JS>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JNS:
            return makeIf<NativeJIT::JccType::JNS>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JP:
            return makeIf<NativeJIT::JccType::JP>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JNP:
            return makeIf<NativeJIT::JccType::JNP>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JL:
            return makeIf<NativeJIT::JccType::JL>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JGE:
            return makeIf<NativeJIT::JccType::JGE>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JLE:
            return makeIf<NativeJIT::JccType::JLE>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        case NativeJIT::JccType::JG:
            return makeIf<NativeJIT::JccType::JG>( expr, cmpLhs, cmpRhs, trueValue, falseValue );
        default:
            throw std::invalid_argument( "Unsupported JccType comparison operator" );
        }
    }

#endif

    template<typename RealType>
    class Parser
    {
    public:
        using Real = RealType;

        //
        // Constructs a parser for the expression text in src.
        // Allocator and FunctionBuffer are constructor parameters
        // for NativeJIT::Function.
        //
        Parser(std::string const & src,
               Allocator& allocator,
               FunctionBuffer& code);

        //
        // Compiles the expression, then invokes the resulting
        // function.
        //
        Real Evaluate( PerfStat &stat );


        //
        // ParseError records the character position and cause of an error
        // during parsing.
        //
        class ParseError : public std::runtime_error
        {
        public:
            ParseError(char const * message, size_t position);

            friend std::ostream& operator<< (std::ostream &out, const ParseError &e);

        private:
            // Character position where error occurred.
            size_t m_position;
        };

    private:
        // Parses an expression of the form
        // EXPRESSION:
        //   SUM
        NativeJIT::Node<RealType>& Parse();

        // Parses expressions of form
        // SUM:
        //   PRODUCT ('+' PRODUCT)*
        //   PRODUCT ('-' PRODUCT)*
        NativeJIT::Node<RealType>& ParseSum();

        // Parses expressions of the form
        // PRODUCT:
        //   TERM ('*' TERM)*
        NativeJIT::Node<RealType>& ParseProduct();

        // Parses expressions of the form
        // TERM:
        //   (SUM)
        //   FLOAT
        //   SYMBOL
        //   "sqrt" '(' SUM ')'
        //   "ifNotZero" '(' SUM ',' SUM ',' SUM  ')'
        //   "if" '(' IF ')'
        NativeJIT::Node<RealType>& ParseTerm();

        // Parses expression of the form
        // IF:
        //   SUM ("<"|"<="|"=="|"!="|">="|">") IFVALUES
        NativeJIT::Node<Real> &ParseIf();

        // Parses expression of the form
        // "<"|"<="|"=="|"!="|">="|">"
        // Returns JccType suitable for floating-point comparison.
        NativeJIT::JccType ParseCmpOp();

        // Parses expressions of the form
        // FLOAT:
        //   [ '+' | '-' ] (DIGIT)* [ '.' DIGIT*] [ ('e' | 'E') [ '+' | '-' ] DIGIT* ]
        Real ParseFloat();

        // Parses expressions of the form
        // SYMBOL: ALPHA (ALPHA | DIGIT)*
        std::string ParseSymbol();

        // Returns true if current position is the first character of a floating
        // point number.
        bool IsFirstCharOfFloat(char c);

        // Advances the current position past whitespace (space, tab, carriage
        // return, newline).
        void SkipWhite();

        // Attempts to advance past next character.
        // Throws a ParseError if the next character is not equal to
        // the expected character.
        void Consume(char expected);

        // Advances past the next character.
        // Returns the character or '\0' if at the end of the stream.
        char GetChar();

        // Returns the next character without advancing.
        // Returns '\0' if at the end of the stream.
        char PeekChar();

        // Source string to be parsed.
        std::string const m_src;

        // Current position of parser in the m_src.
        size_t m_currentPosition;

        // NativeJIT Function used to build and compile parsed expression.
        Function<Real> m_expression;
    };


    template<typename RealType>
    Parser<RealType>::Parser(std::string const & src,
                             Allocator& allocator,
                             FunctionBuffer& code)
        : m_src(src),
          m_currentPosition(0),
          m_expression(allocator, code)
    {
    }

    template<typename Time>
    static long long elapsedNano( Time start, Time end )
    {
        return std::chrono::duration<long long, std::nano>( end - start).count();
    }

    template<typename RealType>
    RealType Parser<RealType>::Evaluate( PerfStat &stat )
    {
        auto startTime = std::chrono::high_resolution_clock::now();
        auto& root = Parse();
        auto parseTime = std::chrono::high_resolution_clock::now();
//        m_expression.EnableDiagnostics( std::cout );
        m_expression.Compile(root);
        auto compileTime = std::chrono::high_resolution_clock::now();
        auto function = m_expression.GetEntryPoint();
        auto result = function();
        auto endTime = std::chrono::high_resolution_clock::now();

        stat.parseNano_ = elapsedNano( startTime, parseTime );
        stat.compileNano_ = elapsedNano( parseTime, compileTime );
        stat.evalNano_ = elapsedNano( compileTime, endTime );

        return result;
    }


    template<typename RealType>
    NativeJIT::Node<RealType>& Parser<RealType>::Parse()
    {
        auto& expression = ParseSum();

        SkipWhite();
        if (PeekChar() != '\0')
        {
            throw ParseError("Syntax error.", m_currentPosition);
        }

        return expression;
    }

    template<typename RealType>
    NativeJIT::Node<RealType>& Parser<RealType>::ParseSum()
    {
        auto& left = ParseProduct();

        SkipWhite();
        if (PeekChar() == '+')
        {
            GetChar();
            auto& right = ParseProduct();

            return m_expression.Add(left, right);
        }
        else if (PeekChar() == '-')
        {
            GetChar();
            auto& right = ParseProduct();

            return m_expression.Sub(left, right);
        }
        else
        {
            return left;
        }
    }


    template<typename RealType>
    NativeJIT::Node<RealType>& Parser<RealType>::ParseProduct()
    {
        auto& left = ParseTerm();

        SkipWhite();
        if (PeekChar() == '*')
        {
            GetChar();
            auto& right = ParseSum();

            return m_expression.Mul(left, right);
        }
        else
        {
            return left;
        }
    }


    template<typename RealType>
    NativeJIT::Node<RealType>& Parser<RealType>::ParseTerm()
    {
        SkipWhite();

        char next = PeekChar();
        if (next == '(')
        {
            GetChar();

            auto& result = ParseSum();

            SkipWhite();
            Consume(')');

            return result;
        }
        else if (IsFirstCharOfFloat(next))
        {
            Real f = ParseFloat();
            return m_expression.Immediate(f);
        }
        else if (isalpha(next))
        {
            std::string symbol = ParseSymbol();
            if (symbol.compare("e") == 0)
            {
                // 'e' denotes Euler's number.
                const Real e = static_cast<Real>(exp(1));
                return m_expression.Immediate(e);
            }
            else if (symbol.compare("pi") == 0)
            {
                // 'pi' denotes the mathematical constant pi.
                const Real pi = static_cast<Real>(atan(1) * 4);
                return m_expression.Immediate(pi);
            }
            else if (symbol.compare( "sqrt" ) == 0)
            {
                Consume( '(' );
                auto& parameter = ParseSum();
                Consume( ')' );
                auto & sqrtFunction = m_expression.Immediate( sqrtf );
                return m_expression.Call( sqrtFunction, parameter );
            }
            else if (symbol.compare( "ifNotZero" ) == 0)
            {
                Consume( '(' );
                auto& conditionValue = ParseSum();
                Consume( ',' );
                auto& trueValue = ParseSum();
                Consume( ',' );
                auto& falseValue = ParseSum();
                Consume( ')' );
                return m_expression.IfNotZero( conditionValue, trueValue, falseValue );
            }
            else if (symbol.compare( "if" ) == 0)
            {
                Consume( '(' );
                auto &ifValue = ParseIf();
                Consume( ')' );
                return ifValue;
            }
            else
            {
                // TODO: REVIEW: Check lifetime of c_str() passed to exception constructor.
                std::stringstream message;
                message << "Unknown identifier \"" << symbol << "\".";
                throw ParseError(message.str().c_str(), m_currentPosition);
            }
        }
        else
        {
            throw ParseError("Expected a number, symbol or parenthesized expression.",
                             m_currentPosition);
        }
    }


    template<typename RealType>
    NativeJIT::Node<RealType>& Parser<RealType>::ParseIf()
    {
        auto& lhs = ParseSum();
        SkipWhite();
        auto cmpOp = ParseCmpOp();
        SkipWhite();
        auto& rhs = ParseSum();
        SkipWhite();
        Consume( ',' );
        auto& ifTrue = ParseSum();
        SkipWhite();
        Consume( ',' );
        auto& ifFalse = ParseSum();
        SkipWhite();
        return If(m_expression, lhs, cmpOp, rhs, ifTrue, ifFalse );
    }


    template<typename RealType>
    NativeJIT::JccType Parser<RealType>::ParseCmpOp()
    {
        // NativeJIT floating point comparison is done using the comiss x86 instruction
        // http://x86.renejeschke.de/html/file_module_x86_id_44.html
        // ZF = 1 => equal
        // ZF = 0 & CF = 1 => less than
        // ZF = 0 & CF = 0 => greater than
        //
        // jb   CF = 1              less than
        // jbe  CF = 1 || ZF = 1    less than or equal
        // je   ZF = 1              equal
        // jbe  ZF = 0              not equal
        // ja   CF = 0 && ZF = 0    greater than
        // jae  CF = 0 || ZF = 1    greater than or equal
        //
        // Not applicable to Real (use flags other ZF/CF set by comiss)
        // jg, jge, jl, jle, jo, js, jns
        // 
        // Visual Studio debugging note: Flags register show with different name.
        // See https://msdn.microsoft.com/en-us/library/kwydd1t7(v=vs.85).aspx
        // PL = SF
        // OF = OV
        // ZF = ZR
        // CY = CF

        auto startPosition = m_currentPosition;
        char first = GetChar();
        switch (first)
        {
        case '<':
            if (PeekChar() == '=')
            {
                GetChar();
                return NativeJIT::JccType::JBE;
            }
            return NativeJIT::JccType::JB;
        case '>':
            if (PeekChar() == '=')
            {
                GetChar();
                return NativeJIT::JccType::JAE;
            }
            return NativeJIT::JccType::JA;
        case '=':
            Consume( '=' );
            return NativeJIT::JccType::JE;
        case '!':
            Consume( '=' );
            return NativeJIT::JccType::JNE;
        };
        throw ParseError( "Expected comparison operator", startPosition );
    }


    template<typename RealType>
    RealType Parser<RealType>::ParseFloat()
    {
        // s will hold a string of floating point number characters that will
        // eventually be passed to stof().
        std::string s;

        SkipWhite();

        //
        // Gather in s the longest possible sequence of characters that will
        // parse as a floating point number.
        //

        auto startPosition = m_currentPosition;

        // Optional leading '+' or '-'.
        if (PeekChar() == '+' || PeekChar() == '-')
        {
            s.push_back(GetChar());
        }

        // Optional mantissa left of '.'
        while (isdigit(PeekChar()))
        {
            s.push_back(GetChar());
        }

        // Optional portion of mantissa right of '.'.
        if (PeekChar() == '.')
        {
            s.push_back(GetChar());
            while (isdigit(PeekChar()))
            {
                s.push_back(GetChar());
            }
        }

        // Optional exponent.
        if (PeekChar() == 'e' || PeekChar() == 'E')
        {
            s.push_back(GetChar());

            // Optional '+' or '-' before exponent.
            if (PeekChar() == '+' || PeekChar() == '-')
            {
                s.push_back(GetChar());
            }

            if (!isdigit(PeekChar()))
            {
                throw ParseError("Expected exponent in floating point constant.",
                                 m_currentPosition);
            }

            while (isdigit(PeekChar()))
            {
                s.push_back(GetChar());
            }
        }

        // Parse s into a floating point value.
        try
        {
            return stof(s);
        }
        catch ( const std::exception& )
        {
            throw ParseError( "Invalid floating point value", startPosition );
        }
    }


    template<typename RealType>
    std::string Parser<RealType>::ParseSymbol()
    {
        std::string symbol;

        SkipWhite();
        if (!isalpha(PeekChar()))
        {
            throw ParseError("Expected alpha character at beginning of symbol.",
                             m_currentPosition);
        }
        while (isalnum(PeekChar()))
        {
            symbol.push_back(GetChar());
        }
        return symbol;
    }


    template<typename RealType>
    bool Parser<RealType>::IsFirstCharOfFloat(char c)
    {
        return isdigit(c) || (c == '-') || (c == '+') || (c == '.');
    }


    template<typename RealType>
    void Parser<RealType>::SkipWhite()
    {
        while (isspace(PeekChar()))
        {
            GetChar();
        }
    }


    template<typename RealType>
    void Parser<RealType>::Consume(char c)
    {
        if (PeekChar() != c)
        {
            // TODO: REVIEW: Check lifetime of c_str() passed to exception constructor.
            std::stringstream message;
            message << "Expected '" << c << "'.";
            throw ParseError(message.str().c_str(), m_currentPosition);
        }
        else
        {
            GetChar();
        }
    }

    
    template<typename RealType>
    char Parser<RealType>::GetChar()
    {
        char result = PeekChar();
        if (result != '\0')
        {
            ++m_currentPosition;
        }
        return result;
    }


    template<typename RealType>
    char Parser<RealType>::PeekChar()
    {
        if (m_currentPosition >= m_src.length())
        {
            return '\0';
        }
        else
        {
            return m_src[m_currentPosition];
        }
    }


    template<typename RealType>
    Parser<RealType>::ParseError::ParseError(char const * message, size_t position)
        : std::runtime_error(message),
          m_position(position)
    {
    }


    template<typename RealType>
    std::ostream& operator<< (std::ostream &out, const typename Parser<RealType>::ParseError &e)
    {
        out << std::string(e.m_position, ' ') << '^' << std::endl;
        out << "Parser error (position = " << e.m_position << "): ";
        out << e.what();
        out << std::endl;
        return out;
    }
}


/******************************************************************************
 *
 * Test() runs a number of test cases for the parser.
 * It prints a summary of each case's input and output and either "OK"
 * or "FAILED" depending on whether the test succeeded or failed.
 *
 * Returns true if all tests pass. Otherwise returns false.
 *
 ******************************************************************************/
bool Test()
{
    class TestCase
    {
    public:
        // TODO: REVIEW: Passing double here instead of float so that I don't
        // have to type the 'f' character at the end of every floating point
        // constant. Consider just changing entire example to use doubles.
        TestCase(char const * input, double output)
            : m_input(input),
            m_output(static_cast<float>(output))
        {
        }


        bool Run(std::ostream& output, Allocator& allocator, FunctionBuffer& code)
        {
            bool succeeded = false;

            output << "\"" << m_input << "\" ==> ";
            Examples::PerfStat perfStat;
            try {
                using Real = float;
                Examples::Parser<Real> parser(m_input, allocator, code);
                Real result = parser.Evaluate( perfStat );

                output << result;

                if (result == m_output)
                {
                    succeeded = true;
                    output << " OK";
                }
                else
                {
                    output << " FAILED: expected " << m_output;
                }
            } catch (...) {
                output << "FAILED: exception.";
            }

            output << std::endl;
            perfStat.print();

            return succeeded;
        }

    private:
        char const * m_input;
        float m_output;
    };


    TestCase cases[] =
    {
        // Constants
        TestCase("1", 1.0),
        TestCase("1.234", 1.234),
        TestCase(".1", 0.1),
        TestCase("-2", -2.0),
        TestCase("-.1", -0.1),
        TestCase("1e9", 1e9),
        TestCase("2e-8", 2e-8),
        TestCase("3e+7", 3e+7),
        TestCase("456.789e+5", 456.789e+5),

        // Symbols
        TestCase("e", static_cast<float>(exp(1))),
        TestCase("pi", static_cast<float>(atan(1) * 4)),

        // Addition
        TestCase("1+2", 3.0),
        TestCase("3+e", 3.0 + static_cast<float>(exp(1))),

        // Subtraction
        TestCase("4-5", -1.0),

        // Multiplication
        TestCase("2*3", 6.0),

        // Parenthesized expressions
        TestCase("(3+4)", 7.0),
        TestCase("(3+4)*(2+3)", 35.0),

        // Combinations
        TestCase("1+-2", -1.0),     // Addition combined with unary negation.

        // White space
        TestCase("\t 1  + ( 2 * 10 )    ", 21.0),

        // sqrt
        TestCase("sqrt(4)", 2.0),
        TestCase("sqrt((3+4)*(2+3))", sqrtf(35)),

        // ifNotZero
        TestCase( "ifNotZero(1, 2, 3)", 2.0 ),
        TestCase( "ifNotZero(0, 2, 3)", 3.0 ),

        // if
        TestCase( "if(1 == 1, 3, 4)", 3.0 ),
        TestCase( "if(1 == 2, 3, 4)", 4.0 ),
        TestCase( "if(1 != 2, 3, 4)", 3.0 ),
        TestCase( "if(1 != 1, 3, 4)", 4.0 ),
        TestCase( "if(1 < 2, 3, 4)", 3.0 ),
        TestCase( "if(5 < 2, 3, 4)", 4.0 ),
        TestCase( "if(2 < 2, 3, 4)", 4.0 ),
        TestCase( "if(1 <= 2, 3, 4)", 3.0 ),
        TestCase( "if(1 <= 1, 3, 4)", 3.0 ),
        TestCase( "if(5 <= 2, 3, 4)", 4.0 ),
        TestCase( "if(5 >= 2, 3, 4)", 3.0 ),
        TestCase( "if(1 >= 2, 3, 4)", 4.0 ),
        TestCase( "if(2 >= 2, 3, 4)", 3.0 ),
        TestCase( "if(5 > 2, 3, 4)", 3.0 ),
        TestCase( "if(1 > 2, 3, 4)", 4.0 ),
        TestCase( "if(2 > 2, 3, 4)", 4.0 ),

        // VS Flags (https://msdn.microsoft.com/en-us/library/kwydd1t7(v=vs.85).aspx)
        // PL = SF
        // OF = OV
        // ZF = ZR
        // CY = CF
        // jg: PL = OV and ZR = 0
        // 1 < 2: OV = 0 UP = 0 EI = 1 PL = 0 ZR = 0 AC = 0 PE = 0 CY = 1
        // 5 < 2: OV = 0 UP = 0 EI = 1 PL = 0 ZR = 0 AC = 0 PE = 0 CY = 0 
        // => only difference in on CY, the carry flag
    };


    ExecutionBuffer codeAllocator(8192);
    Allocator allocator(8192);
    FunctionBuffer code(codeAllocator, 8192);
    //code.EnableDiagnostics( std::cout );

    bool success = true;
    for (size_t i = 0; i < sizeof(cases) / sizeof(TestCase); ++i)
    {
        allocator.Reset();
        codeAllocator.Reset();
//        if ( i == (sizeof( cases ) / sizeof( TestCase )-1) )
        {
            success &= cases[i].Run(std::cout, allocator, code);
        }
    }

    return success;
}


int main( int argc, const char *argv[] )
{
#ifndef NATIVEJIT_WITH_AFL
    std::cout << "Running test cases ..." << std::endl;
    bool success = Test();
    if (success)
    {
        std::cout << "All tests succeeded." << std::endl;
    }
    else
    {
        std::cout << "One or more tests failed." << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Type an expression and press return to evaluate." << std::endl;
    std::cout << "Enter an empty line to exit." << std::endl;
#endif

    ExecutionBuffer codeAllocator(8192);
    Allocator allocator(8192);
    FunctionBuffer code(codeAllocator, 8192);
    //code.EnableDiagnostics( std::cout ); // uncomment to see assembly generated for the function
    std::string prompt(">> ");

    std::istream *inputStream = &std::cin;
    std::ifstream ifs;
    if ( argc > 1 )
    {
        ifs.open( argv[1] );
        if ( !ifs.good() )
        {
            std::cout << "Failed to open file '" << argv[1] << "'." << std::endl;
            return 1;
        }
        inputStream = &ifs;
    }

#ifdef NATIVEJIT_WITH_AFL
    while ( __AFL_LOOP( 1000 ) )
#endif
    for (;;)
    {
        allocator.Reset();
        codeAllocator.Reset();

        std::string line;
        std::cout << prompt << std::flush;
        std::getline(*inputStream, line);

        // TODO: Should really see if line is completely blank.
        // Blank lines cause the parser to crash.
        if (line.length() == 0)
        {
            break;
        }

        using Real = float;
        try
        {
            Examples::Parser<Real> parser(line, allocator, code);
            Examples::PerfStat perfStat;
            float result = parser.Evaluate( perfStat );
            std::cout << result << std::endl;
#ifndef NATIVEJIT_WITH_AFL
            perfStat.print();
#endif
        }
        catch (Examples::Parser<Real>::ParseError& e)
        {
            std::cout << std::string(prompt.length(), ' ');
            std::cout << e;
        }
    }

    return 0;
}
