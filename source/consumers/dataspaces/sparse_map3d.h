#pragma once

#include "dataspace.h"

namespace DAQuiri
{

class SparseMap3D : public Dataspace
{
  public:
    SparseMap3D();
    SparseMap3D* clone() const override
    { return new SparseMap3D(*this); }

    void clear() override;
    void add(const Entry&) override;
    void add_one(const Coords&) override;
    PreciseFloat get(const Coords&) const override;
    EntryList range(std::vector<Pair> list) const override;
    void recalc_axes() override;

    void save(std::ostream& os) override;

  protected:
    using tripple = std::tuple<uint16_t,uint16_t,uint16_t>;
    using SpectrumMap3D = std::map<tripple, PreciseFloat>;

    //the data itself
    SpectrumMap3D spectrum_;
    PreciseFloat total_count_ {0};
    uint16_t max0_ {0};
    uint16_t max1_ {0};
    uint16_t max2_ {0};

    inline void bin_pair(const uint16_t& x, const uint16_t& y, const uint16_t& z,
                         const PreciseFloat& count)
    {
      spectrum_[tripple(x,y,z)] += count;
      total_count_ += count;
      max0_ = std::max(max0_, x);
      max1_ = std::max(max1_, y);
      max2_ = std::max(max2_, z);
    }

    inline void bin_one(const uint16_t& x, const uint16_t& y, const uint16_t& z)
    {
      spectrum_[tripple(x,y,z)] ++;
      total_count_ ++;
      max0_ = std::max(max0_, x);
      max1_ = std::max(max1_, y);
      max2_ = std::max(max2_, z);
    }

    bool is_symmetric();

    void fill_list(EntryList &result,
                   size_t min0, size_t max0,
                   size_t min1, size_t max1,
                   size_t min2, size_t max2) const;

    void data_save(hdf5::node::Group) const override;
    void data_load(hdf5::node::Group) override;

    std::string data_debug(const std::string& prepend) const override;
};

}
