///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2014) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * FeudalObjectManager.h
 *
 *  Created on: Sep 16, 2014
 *      Author: tsharpe
 */
#ifndef FEUDAL_OBJECTMANAGER_H_
#define FEUDAL_OBJECTMANAGER_H_

#include "feudal/BinaryStream.h"
#include "feudal/MasterVec.h"
#include "system/Assert.h"
#include "system/file/File.h"

template <class T>
class ObjectManagerImpl
{
protected:
    void loadImpl( File const& file, T* pObj )
    { BinaryReader::readFile(file,pObj); }
    void storeImpl( File const& file, T const& obj )
    { BinaryWriter::writeFile(file,obj); }
};

template <class T>
class ObjectManagerImpl<MasterVec<T>>
{
protected:
    void loadImpl( File const& file, MasterVec<T>* pObj )
    { pObj->ReadAll(String(file.c_str())); }
    void storeImpl( File const& file, MasterVec<T> const& obj )
    { obj.WriteAll(String(file.c_str())); }
};

template <class T>
class ObjectManager : ObjectManagerImpl<T>
{
public:
    ObjectManager( String const& fileName ) : mFile(fileName), mpObj(nullptr) {}
    ObjectManager( ObjectManager const& )=delete;
    ObjectManager& operator=( ObjectManager const& )=delete;
    ~ObjectManager() { delete mpObj; }

    T& create() { delete mpObj; mpObj = new T; return *mpObj; }

    T const& load()
    { if ( !mpObj ) { mpObj=new T; this->loadImpl(mFile,mpObj); }
      return *mpObj; }

    void unload() { delete mpObj; mpObj = nullptr; }
    void store() { ForceAssert(mpObj); this->storeImpl(mFile,*mpObj); }
    bool onDisk() { return mFile.exists(true); }
    void remove() { mFile.remove(); }

private:
    File mFile;
    T* mpObj;
};

#endif /* FEUDAL_OBJECTMANAGER_H_ */
