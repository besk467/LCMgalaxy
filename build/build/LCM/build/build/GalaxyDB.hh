/*
 GalaxyDB
 DB of galaxy file information

 Created by: Noah Oblath
 Date: Jan 8, 2013
 */

#include <map>
#include <string>


class GalaxyDB
{
public:
	typedef std::map< std::string, std::string > GalaxyFiles;
private:
	typedef std::map< std::string, GalaxyFiles > FileMap;

public:
	GalaxyDB();
	GalaxyDB(const std::string& filename, const std::string& ltName="Label", const std::string& dirName="Dir", const std::string& dataPrefix="Data");
	~GalaxyDB();

	bool LoadFile(const std::string& filename);

	bool DBIsReady() const;

	bool HasGalaxy(const std::string& galaxyLabel) const;
    bool HasGalaxyData(const std::string& galaxyDataLabel) const;

	const std::string& GetFile(const std::string& galaxyLabel, const std::string& fileType) const;
	const GalaxyFiles& GetFiles(const std::string& galaxyLabel) const;

	const std::string& GetDataFile(const std::string& galaxyLabel, const std::string& fileType) const;
	const GalaxyFiles& GetDataFiles(const std::string& galaxyLabel) const;
    
	const std::string& FileNotFound() const;

	const std::string& LabelTokenName() const;
	void SetLabelTokenName(const std::string& name);

	const std::string& DirTokenName() const;
	void SetDirTokenName(const std::string& name);

private:
	bool AddGalaxy(const std::string& line);

	std::string Trim(const std::string& line, const std::string& whitespace = " \t");

	FileMap galaxyFiles;
    FileMap galaxyDataFiles;

	std::string fileNotFound;

	std::string labelTokenName;
	std::string dirTokenName;
    std::string dataTokenPrefix;
};
