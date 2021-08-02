// an attempt to parse the entire database into C++ "Block" objects
// tedious and incomplete
// Author: Mason Housenga

#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <ctime>
#include <algorithm>

const std::vector<int> blockTypes{123, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24};
const std::vector<std::string> statuses{"0", "1", "2", "3", "4", "5", "5a", "5b", "5c", "5m", "5r", "6", "7", "8"};

class Block {
  public:
  Block(std::vector<std::string> row_data, int from_sheet, int from_row) : 
    dbn(row_data[0]), block(row_data[2]), status(row_data[3]), shipment(row_data[4]), shipmentDate(row_data[5]),
    comment(row_data[6]), powder(row_data[14]), bucket(row_data[15]), moldSeries(row_data[16]), emptyMoldMass(row_data[17]),
    filledMoldMass(row_data[18]), tungstenMass(row_data[19]), fiberLoader(row_data[20]), tungstenFillingDate(row_data[21]), tungstenFiller(row_data[22]),
    epoxyBatch(row_data[23]), resinMass(row_data[24]), hardenerMass(row_data[25]), pottingNotes(row_data[26]), epoxyFillingTime(row_data[27]),
    epoxyFillingDate(row_data[28]), epoxyPreparer(row_data[29]), machiningDate(row_data[31]), machinist(row_data[32]), l(row_data[33]),
    bt(row_data[34]), bb(row_data[35]), bh(row_data[36]), st(row_data[37]), sb(row_data[38]),
    sh(row_data[39]), volume(row_data[40]), dimensionTester(row_data[41]), mass(row_data[42]), massTester(row_data[43]),
    density(row_data[44]), densityDate(row_data[45]), goodEnd(row_data[47]), fiberPercentage(row_data[49]), tower1(row_data[55]),
    tower2(row_data[56]), tower3(row_data[57]), tower4(row_data[58]), missingRow(row_data[59]), thirteenHoles(row_data[60]),
    lightTransTester(row_data[63]), lightTransDate(row_data[64]), scintillation(row_data[65]), scintRatio(row_data[67]), scintTester(row_data[68]),
    scintDate(row_data[69]), natLightDate(row_data[70]), natLightTester(row_data[71]), densityCheckedBy(row_data[73]), lightTransCheckedBy(row_data[74]),
    natLightCheckedBy(row_data[75]), pregrade(row_data[80]) {
    sheet = from_sheet;
    row = from_row;
  }
  // custom member classes for all parameters
  class BlockValue {
    public:
      BlockValue(std::string data = "") {
        hasData = dataPresent(data);
        str = data;
      }
      std::string str;
      bool hasData;
  };

  class BlockInt : public BlockValue {
    public:
      BlockInt(std::string data = "") : BlockValue(data) {
        try {
          int value = std::stoi(data);
          valid = true;
        } catch(std::exception& e) {
          // failed to convert to int
          valid = false;
        }
      }
      bool valid;
      int value;
  };

  class BlockDouble : public BlockValue {
    public:
      BlockDouble(std::string data = "") : BlockValue(data) {
        try {
          double value = std::stod(data);
          valid = true;
        } catch (std::exception& e) {
          // failed to convert to double
          valid = false;
        }
      }
      bool valid;
      double value;
  };

  class BlockUSDate : public BlockValue {
    public:
      BlockUSDate(std::string data = "") : BlockValue(data) {
        std::vector<std::string> split = split_string(data, "/");
        if (split.size() == 3) {
          try {
            int m = std::stoi(split[0]);
            int d = std::stoi(split[1]);
            int y = std::stoi(split[2]);
            struct std::tm date = {0,0,0,d,m-1,y-1900};
            time = std::mktime(&date);
            if (time != (std::time_t)(-1)) {
              valid = true;
            } else {
              valid = false;
            }
          } catch (std::exception& e) {
            // failed to convert to date
            valid = false;
          }
        } else {
          valid = false;
        }
      }
      bool valid;
      std::time_t time;
  };

  class BlockNumDate : public BlockValue {
    public:
      BlockNumDate(std::string data = "") : BlockValue(data) {
        std::vector<std::string> split = split_string(data, ".");
        if (split.size() == 1) {
          split.push_back("0");
        }
        if (split.size() == 2) {
          try {
            std::string prefix = split[0];
            instance = std::stoi(split[1]);
            if (prefix.length() == 8) {
              int y = std::stoi(prefix.substr(0, 4));
              int m = std::stoi(prefix.substr(4, 2));
              int d = std::stoi(prefix.substr(6, 2));
              struct std::tm date = {0,0,0,d,m-1,y-1900};
              time = std::mktime(&date);
              if (time != (std::time_t)(-1)) {
                double value = std::stod(data);
                valid = true;
              } else {
                valid = false;
              }
            } else {
              valid = false;
            }
          } catch (std::exception& e) {
            // failed to convert to double
            valid = false;
          }
        } else {
          valid = false;
        }
      }
      bool valid;
      double num;
      double instance;
      std::time_t time;
  };
  class dbn : public BlockValue {
    public:
      dbn(std::string data = "") : BlockValue(data) {}
      bool has = hasData;
      std::string get() {
        return str;
      }
    private:
  };
  class block : public BlockInt {
    public:
      block(std::string data = "") : BlockInt(data) {
        has = valid && std::find(blockTypes.begin(), blockTypes.end(), value) != blockTypes.end();
      }
      bool has;
      int get() {
        return value;
      }
    private:
  };
  class status : public BlockValue {
    public:
      status(std::string data = "") : BlockValue(data) {
        has = std::find(statuses.begin(), statuses.end(), str) != statuses.end();
      }
      bool has;
      std::string get() {
        return str;
      }
    private:
  };
  class shipment : public BlockInt {
    public:
      shipment(std::string data = "") : BlockInt(data) {
        has = valid;
      }
      bool has;
      int get() {
        return value;
      }
    private:
  };
  class shipmentDate : public BlockUSDate {
    public:
      shipmentDate(std::string data = "") : BlockUSDate(data) {
        has = valid;
      }
      bool has;
      std::string get() {
        return str;
      }
    private:
  };
  class comment : public BlockValue {
    public:
      comment(std::string data = "") : BlockValue(data) {
        has = hasData;
      }
      bool has;
      std::string get() {
        return str;
      }
    private:
  };
  class powder : public BlockValue {
    public:
      powder(std::string data = "") : BlockValue(data) {
        has = hasData;
      }
      bool has;
      std::string get() {
        return str;
      }
    private:
  };
  class bucket : public BlockValue {
    public:
      bucket(std::string data = "") : BlockValue(data) {
        has = hasData;
      }
      bool has;
      std::string get() {
        return str;
      }
    private:
  };
  class moldSeries : public BlockValue {
    public:
      moldSeries(std::string data = "") : BlockValue(data) {
        has = hasData;
      }
      bool has;
      std::string get() {
        return str;
      }
    private:
  };
  class emptyMoldMass : public BlockDouble {
    public:
      emptyMoldMass(std::string data = "") : BlockDouble(data) {
        has = valid;
      }
      bool has;
      int get() {
        return value;
      }
    private:
  };
  class filledMoldMass : public BlockDouble {
    public:
      filledMoldMass(std::string data = "") : BlockDouble(data) {
        has = valid;
      }
      bool has;
      double get() {
        return value;
      }
    private:
  };
  class tungstenMass : public BlockDouble {
    public:
      tungstenMass(std::string data = "") : BlockDouble(data) {
        has = valid;
      }
      bool has;
      double get() {
        return value;
      }
    private:
  };
  class fiberLoader : public BlockValue {
    public:
      fiberLoader(std::string data = "") : BlockValue(data) {}
    private:
  };
  class tungstenFillingDate : public BlockUSDate {
    public:
      tungstenFillingDate(std::string data = "") : BlockUSDate(data) {}
    private:
  };
  class tungstenFiller : public BlockValue {
    public:
      tungstenFiller(std::string data = "") : BlockValue(data) {}
    private:
  };
  class epoxyBatch : public BlockInt {
    public:
      epoxyBatch(std::string data = "") : BlockInt(data) {}
    private:
  };
  class resinMass : public BlockDouble {
    public:
      resinMass(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class hardenerMass : public BlockDouble {
    public:
      hardenerMass(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class pottingNotes : public BlockValue {
    public:
      pottingNotes(std::string data = "") : BlockValue(data) {}
    private:
  };
  class epoxyFillingTime : public BlockDouble {
    public:
      epoxyFillingTime(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class epoxyFillingDate : public BlockUSDate {
    public:
      epoxyFillingDate(std::string data = "") : BlockUSDate(data) {}
    private:
  };
  class epoxyPreparer : public BlockValue {
    public:
      epoxyPreparer(std::string data = "") : BlockValue(data) {}
    private:
  };
  class machiningDate : public BlockUSDate {
    public:
      machiningDate(std::string data = "") : BlockUSDate(data) {}
    private:
  };
  class machinist : public BlockValue {
    public:
      machinist(std::string data = "") : BlockValue(data) {}
    private:
  };
  class l : public BlockDouble {
    public:
      l(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class bt : public BlockDouble {
    public:
      bt(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class bb : public BlockDouble {
    public:
      bb(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class bh : public BlockDouble {
    public:
      bh(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class st : public BlockDouble {
    public:
      st(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class sb : public BlockDouble {
    public:
      sb(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class sh : public BlockDouble {
    public:
      sh(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class volume : public BlockDouble {
    public:
      volume(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class dimensionTester : public BlockValue {
    public:
      dimensionTester(std::string data = "") : BlockValue(data) {}
    private:
  };
  class mass : public BlockDouble {
    public:
      mass(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class massTester : public BlockValue {
    public:
      massTester(std::string data = "") : BlockValue(data) {}
    private:
  };
  class density : public BlockDouble {
    public:
      density(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class densityDate : public BlockNumDate {
    public:
      densityDate(std::string data = "") : BlockNumDate(data) {}
    private:
  };
  class goodEnd : public BlockValue {
    public:
      goodEnd(std::string data = "") : BlockValue(data) {}
    private:
  };
  class fiberPercentage : public BlockDouble {
    public:
      fiberPercentage(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class tower1 : public BlockDouble {
    public:
      tower1(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class tower2 : public BlockDouble {
    public:
      tower2(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class tower3 : public BlockDouble {
    public:
      tower3(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class tower4 : public BlockDouble {
    public:
      tower4(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class missingRow : public BlockValue {
    public:
      missingRow(std::string data = "") : BlockValue(data) {}
    private:
  };
  class thirteenHoles : public BlockValue {
    public:
      thirteenHoles(std::string data = "") : BlockValue(data) {}
    private:
  };
  class lightTransTester : public BlockValue {
    public:
      lightTransTester(std::string data = "") : BlockValue(data) {}
    private:
  };
  class lightTransDate : public BlockNumDate {
    public:
      lightTransDate(std::string data = "") : BlockNumDate(data) {}
    private:
  };
  class scintillation : public BlockDouble {
    public:
      scintillation(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class scintRatio : public BlockDouble {
    public:
      scintRatio(std::string data = "") : BlockDouble(data) {}
    private:
  };
  class scintTester : public BlockValue {
    public:
      scintTester(std::string data = "") : BlockValue(data) {}
    private:
  };
  class scintDate : public BlockNumDate {
    public:
      scintDate(std::string data = "") : BlockNumDate(data) {}
    private:
  };
  class natLightDate : public BlockNumDate {
    public:
      natLightDate(std::string data = "") : BlockNumDate(data) {}
    private:
  };
  class natLightTester : public BlockValue {
    public:
      natLightTester(std::string data = "") : BlockValue(data) {}
    private:
  };
  class densityCheckedBy : public BlockValue {
    public:
      densityCheckedBy(std::string data = "") : BlockValue(data) {}
    private:
  };
  class lightTransCheckedBy : public BlockValue {
    public:
      lightTransCheckedBy(std::string data = "") : BlockValue(data) {}
    private:
  };
  class natLightCheckedBy : public BlockValue {
    public:
      natLightCheckedBy(std::string data = "") : BlockValue(data) {}
    private:
  };
  class pregrade : public BlockValue {
    public:
      pregrade(std::string data = "") : BlockValue(data) {}
    private:
  };
  int sheet;
  int row;
  Block::dbn dbn;
  Block::block block;
  Block::status status;
  Block::shipment shipment;
  Block::shipmentDate shipmentDate;
  Block::comment comment;
  Block::powder powder;
  Block::bucket bucket;
  Block::moldSeries moldSeries;
  Block::emptyMoldMass emptyMoldMass;
  Block::filledMoldMass filledMoldMass;
  Block::tungstenMass tungstenMass;
  Block::fiberLoader fiberLoader;
  Block::tungstenFillingDate tungstenFillingDate;
  Block::tungstenFiller tungstenFiller;
  Block::epoxyBatch epoxyBatch;
  Block::resinMass resinMass;
  Block::hardenerMass hardenerMass;
  Block::pottingNotes pottingNotes;
  Block::epoxyFillingTime epoxyFillingTime;
  Block::epoxyFillingDate epoxyFillingDate;
  Block::epoxyPreparer epoxyPreparer;
  Block::machiningDate machiningDate;
  Block::machinist machinist;
  Block::l l;
  Block::bt bt;
  Block::bb bb;
  Block::bh bh;
  Block::st st;
  Block::sb sb;
  Block::sh sh;
  Block::volume volume;
  Block::dimensionTester dimensionTester;
  Block::mass mass;
  Block::massTester massTester;
  Block::density density;
  Block::densityDate densityDate;
  Block::goodEnd goodEnd;
  Block::fiberPercentage fiberPercentage;
  Block::tower1 tower1;
  Block::tower2 tower2;
  Block::tower3 tower3;
  Block::tower4 tower4;
  Block::missingRow missingRow;
  Block::thirteenHoles thirteenHoles;
  Block::lightTransTester lightTransTester;
  Block::lightTransDate lightTransDate;
  Block::scintillation scintillation;
  Block::scintRatio scintRatio;
  Block::scintTester scintTester;
  Block::scintDate scintDate;
  Block::natLightDate natLightDate;
  Block::natLightTester natLightTester;
  Block::densityCheckedBy densityCheckedBy;
  Block::lightTransCheckedBy lightTransCheckedBy;
  Block::natLightCheckedBy natLightCheckedBy;
  Block::pregrade pregrade;
  private:
};

std::vector<std::string> split_string(std::string s, std::string delimiter) {
  size_t pos_start = 0, pos_end, delim_len = delimiter.length();
  std::string token;
  std::vector<std::string> res;
  while ((pos_end = s.find (delimiter, pos_start)) != std::string::npos) {
      token = s.substr (pos_start, pos_end - pos_start);
      pos_start = pos_end + delim_len;
      res.push_back (token);
  }
  res.push_back (s.substr (pos_start));
  return res;
}

bool dataPresent(std::string data) {
  return data != "" && data != "#NUM!" && data != "#DIV/0!" && data != "NaN";
}