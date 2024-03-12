#ifndef PARAMETERREADER_HPP_
#define PARAMETERREADER_HPP_

#include "Exception.hpp"
#include "FileFinder.hpp"
#include <unordered_map>
#include <vector>

#define READ_PARAMETER(READER, TYPE, NAME) TYPE const NAME = static_cast<TYPE>(READER.get(#NAME));

#define PARAMETERIZED_TEST(FUNC, DIR, PARAM_FILE) \
    void Test##PARAM_FILE()                       \
    {                                             \
        FUNC(DIR #PARAM_FILE);                    \
    }

inline bool isspace(char charc)
{
    return charc == ' ' || charc == '\t' || charc == '\n' || charc == '\v' || charc == '\f' || charc == '\r';
}

inline bool isdigit(char charc) { return ('0' <= charc && charc <= '9') || charc == '.'; }

inline bool isalpha(char charc)
{
    return ('a' <= charc && charc <= 'z') || ('A' <= charc && charc <= 'Z');
}

inline bool isalnum(char charc) { return isdigit(charc) || isalpha(charc) || charc == '_'; }

class ParameterReader
{
public:
    void parse(std::string const& path)
    {
        m_slices.clear();
        std::ifstream fstream{ path };
        if (!fstream.is_open())
        {
            EXCEPTION("Can't open file " + path);
        }
        char c = 0;
        std::string buf{};
        while (fstream >> std::noskipws >> c)
        {
            if (isalnum(c))
            {
                buf.push_back(c);
            }
            else
            {
                if (!buf.empty())
                {
                    m_slices.push_back(buf);
                    buf.clear();
                }
            }
        }
        fstream.close();
        if (m_slices.size() % 2 != 0)
        {
            EXCEPTION("Configuration file should be key value pairs.");
        }
        try
        {
            for (size_t i = 0; i < m_slices.size(); i += 2)
            {
                double const value = std::stod(m_slices[i + 1]);
                m_data[m_slices[i]] = value;
            }
        }
        catch (std::exception& exp)
        {
            EXCEPTION(exp.what());
        }
    }

    double get(const char* key) const
    {
        try
        {
            return m_data.find(std::string{ key })->second;
        }
        catch (std::exception& exp)
        {
            EXCEPTION(exp.what());
        }
    }

private:
    std::vector<std::string> m_slices{};

    std::unordered_map<std::string, double> m_data{};
};

#endif /*PARAMETERREADER_HPP_*/
