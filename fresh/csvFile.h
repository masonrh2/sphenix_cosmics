#pragma once
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

class csvfile
{
    std::ofstream fs_;
    bool is_first_;
    const std::string separator_;
public:
    csvfile(const std::string filename, const std::string separator = ",")
        : fs_()
        , is_first_(true)
        , separator_(separator)
    {
        fs_.exceptions(std::ios::failbit | std::ios::badbit);
        fs_.open(filename);
    }

    ~csvfile()
    {
        _flush();
        fs_.close();
    }

    void _flush()
    {
        fs_.flush();
    }

    void _endrow()
    {
        fs_ << std::endl;
        is_first_ = true;
    }

    csvfile& operator << (csvfile& (* val)(csvfile&))
    {
        return val(*this);
    }

    csvfile& operator << (const char * val)
    {
        return write(val);
    }

    csvfile& operator << (const std::string & val)
    {
        return write(val);
    }

    template<typename T>
    csvfile& operator << (const T& val)
    {
        return write(val);
    }

    inline static csvfile& endrow(csvfile& file)
    {
        file._endrow();
        return file;
    }

    inline static csvfile& flush(csvfile& file)
    {
        file._flush();
        return file;
    }
private:
    template<typename T>
    csvfile& write (const T& val)
    {
        if (!is_first_)
        {
            fs_ << separator_;
        }
        else
        {
            is_first_ = false;
        }
        fs_ << val;
        return *this;
    }

};