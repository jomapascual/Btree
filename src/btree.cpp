/***
 * @author See Contributors.txt for code contributors and overview of BadgerDB.
 *
 * @section LICENSE
 * Copyright (c) 2012 Database Group, Computer Sciences Department, University of Wisconsin-Madison.
 */

#include "btree.h"
#include "filescan.h"
#include "exceptions/bad_index_info_exception.h"
#include "exceptions/bad_opcodes_exception.h"
#include "exceptions/bad_scanrange_exception.h"
#include "exceptions/no_such_key_found_exception.h"
#include "exceptions/scan_not_initialized_exception.h"
#include "exceptions/index_scan_completed_exception.h"
#include "exceptions/file_not_found_exception.h"
#include "exceptions/end_of_file_exception.h"


//#define DEBUG

namespace badgerdb
{

// -----------------------------------------------------------------------------
// BTreeIndex::BTreeIndex -- Constructor
// -----------------------------------------------------------------------------

BTreeIndex::BTreeIndex(const std::string & relationName,
		std::string & outIndexName,
		BufMgr *bufMgrIn,
		const int attrByteOffset,
		const Datatype attrType)
{
    // Add your code below. Please do not remove this line.
	bufMgr = bufMgrIn;
	//this.attrByteOffset -> attrByteOffset;
	attributeType = attrType;

	std::ostringstream idxStr;
	idxStr << relationName << "." << attrByteOffset;
	std::string indexName = idxStr.str();

	Page *headerPage, *rootPage, *leafPage;
	IndexMetaInfo* meta;
	try {
		file = new BlobFile(indexName, false); // Checks if the file exists

		headerPageNum = file -> getFirstPageNo();
		bufMgr -> readPage(file, headerPageNum, headerPage); // Reads the first page (meta)
		rootPageNum = ((IndexMetaInfo*) headerPage) -> rootPageNo; // Held in the meta page

		bufMgr -> unPinPage(file, headerPageNum, false); 
		bufMgr -> unPinPage(file, rootPageNum, false);
	}
	catch (FileNotFoundException e) {
		bufMgr -> allocPage(file, headerPageNum, headerPage); // Allocate header page
    	bufMgr -> allocPage(file, rootPageNum, rootPage); // Allocate root page

		meta = (IndexMetaInfo*) headerPage; // Fills in meta info
		strcpy(meta -> relationName, relationName.c_str());
		meta -> attrByteOffset = attrByteOffset;
		meta -> attrType = attrType;
		meta -> rootPageNo = rootPageNum;
		
		LeafNodeInt *root = (LeafNodeInt *)rootPage; // Initializes the root
    	root -> rightSibPageNo = 0;

		bufMgr -> unPinPage(file, headerPageNum, false); 
		bufMgr -> unPinPage(file, rootPageNum, false); 

		// Scans Records
		FileScan scan = FileScan(relationName, bufMgr);
		scanExecuting = true;
		RecordId rid;
		try {
			while(1) {
				scan.scanNext(rid);
				std::string fileRecord = scan.getRecord();
				const char *record = fileRecord.c_str();
				insertEntry(record + attrByteOffset, rid);
			}
		}
		catch (EndOfFileException e) {
			bufMgr -> flushFile(file);
		}
	}
}


// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------

BTreeIndex::~BTreeIndex()
{
    // Add your code below. Please do not remove this line.
	scanExecuting = false;
	bufMgr -> flushFile(file);
	delete file;
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

 /**
	 * Insert a new entry using the pair <value,rid>. 
	 * Start from root to recursively find out the leaf to insert the entry in. The insertion may cause splitting of leaf node.
	 * This splitting will require addition of new leaf page number entry into the parent non-leaf, which may in-turn get split.
	 * This may continue all the way upto the root causing the root to get split. If root gets split, metapage needs to be changed accordingly.
	 * Make sure to unpin pages as soon as you can.
   * @param key			Key to insert, pointer to integer/double/char string
   * @param rid			Record ID of a record whose entry is getting inserted into the index.
	**/
void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
    // Add your code below. Please do not remove this line.

	// Create pair
	RIDKeyPair<const void*> newEntry;
	newEntry.set(rid, key);

}

void BTreeIndex::insertRoot(){
	
}

void BTreeIndex::insertLeaf(LeafNodeInt *leafNode, RIDKeyPair<int> ridKey) // Should use parameters like this
{
	if (leafNode -> ridArray[0].page_number == 0) { // Page is empty
		leafNode -> keyArray[0] = ridKey.key;
		leafNode -> ridArray[0] = ridKey.rid;    
	}
	else {
		int i = INTARRAYLEAFSIZE; // Number of keys in the leaf node // Maybe change to leafOccupancy ?
		while((leafNode -> ridArray[i-1].page_number == 0) && i > 0) { // Gets to the end of leafNode
			i--;
		}

		while((leafNode -> keyArray[i-1] > ridKey.key) && i > 0) { // Shifts the previous ridKey pairs
			leafNode -> keyArray[i] = leafNode -> keyArray[i-1];
			leafNode -> ridArray[i] = leafNode -> ridArray[i-1];
			i--;
		}

		leafNode -> keyArray[i] = ridKey.key; // Inserts ridKey to the leafNode
		leafNode -> ridArray[i] = ridKey.rid;
	}
}

void BTreeIndex::insertNonLeaf(NonLeafNodeInt *node, PageKeyPair<int> keyPage)
{
	int i = INTARRAYNONLEAFSIZE; // Maybe change to nodeOccupancy ?
	while((node -> pageNoArray[i] == 0) && i > 0) {
		i--;
	}
	while((node -> keyArray[i-1] > keyPage.key) && i > 0) {
		node -> pageNoArray[i+1] = node -> pageNoArray[i]; // TODO not 100% sure if this is the right way to shift pageNo
		node -> keyArray[i] = node -> pageNoArray[i-1]; // Shifts the previous keyPage pairs
		i--;
	}

	node -> pageNoArray[i+1] = keyPage.pageNo; // Inserts keyPage pair to the node
	node -> keyArray[i] = keyPage.key;
}

  PageKeyPair<int> BTreeIndex::splitLeafNode(LeafNodeInt *oldLeafNode, PageId currPageId, RIDKeyPair<int> ridPair)
  {
	// allocate new leaf sibling
	Page* sibling;
	PageId siblingId; // new leaf sibling page number
	bufMgr -> allocPage(file, siblingId, sibling);
	// convert to proper structure
	LeafNodeInt* siblingNode = (LeafNodeInt*) sibling;

	// current leaf node gets right sibling page number added to it
	if (oldLeafNode -> rightSibPageNo != 0) {
		siblingNode -> rightSibPageNo = oldLeafNode -> rightSibPageNo;
	}
	oldLeafNode -> rightSibPageNo = siblingId;

	// split current leaf into two seperate leaves
	for (int i = 0; i < INTARRAYLEAFSIZE / 2; ++i) {
		// split somehow:
		siblingNode -> keyArray[i] = oldLeafNode -> keyArray[(INTARRAYLEAFSIZE/2) + i]; // Right half of old node gets put into new node
		siblingNode -> ridArray[i] = oldLeafNode -> ridArray[(INTARRAYLEAFSIZE/2) + i];
		oldLeafNode -> keyArray[(INTARRAYLEAFSIZE/2) + i] = 0; // Right half of old node set to 0
		oldLeafNode -> ridArray[(INTARRAYLEAFSIZE/2) + i].page_number = 0;
	}

	// insert pair into the split leaves
	if (ridPair.key < siblingNode -> keyArray[0]) {
		insertLeaf(oldLeafNode, ridPair);
	}
	else {
		insertLeaf(siblingNode, ridPair);
	}

	// calculate new middle key pair:
	PageKeyPair<int> middleKeyPair;
	middleKeyPair.set(siblingId, siblingNode -> keyArray[0]);

	// return recursive to move up tree
	return middleKeyPair;
  }

//splitNonLeaf
//splitRoot

void BTreeIndex::splitNonLeafNode(NonLeafNodeInt *oldNonLeafNode, PageId currPageId, PageKeyPair<int>*&newPair)
 {
	 // allocate new non-leaf sibling
	 Page* sibling;
	 PageId siblingId;
	 bufMgr -> allocPage(file, siblingId, sibling);
	 // convert to proper structure
	 NonLeafNodeInt* siblingNode = (NonLeafNodeInt*) sibling;

	 // middle range
	 int middle = nodeOccupancy / 2;
	 int pushUp = middle;
	 PageKeyPair<int> nodeToPush;

	// check if full, to push up
	 if(nodeOccupancy %2 == 0)
	 {
		 if(newPair -> key < oldNonLeafNode -> keyArray[middle])
		 {
			 pushUp = middle - 1;
		 } else{
			 pushUp = middle;
		 }
	 }
	 nodeToPush.set(siblingId, oldNonLeafNode->keyArray[pushUp]);

	middle = pushUp + 1;
	for(int i = middle; i < nodeOccupancy; i++)
	{
		siblingNode -> keyArray[i - middle] = oldNonLeafNode -> keyArray[i];
		siblingNode -> pageNoArray[i - middle] = oldNonLeafNode -> pageNoArray[i + 1];
		oldNonLeafNode -> pageNoArray[i + 1] = (PageId) 0;
		oldNonLeafNode -> keyArray[i + 1] = 0;
	}

	siblingNode -> level = oldNonLeafNode -> level;
	oldNonLeafNode -> keyArray[pushUp] = 0;
	oldNonLeafNode -> pageNoArray[pushUp] = (PageId) 0;

	if(newPair -> key < siblingNode -> keyArray[0])
	{
		insertNonLeaf(oldNonLeafNode, *newPair);
	} else{
		insertNonLeaf(siblingNode, *newPair);
	}

	newPair = &nodeToPush;
	bufMgr -> unPinPage(file, currPageId, true);
	bufMgr -> unPinPage(file, siblingId, true);

	if(currPageId == rootPageNum)
	{
		updateRootNode(*newPair, currPageId);
	}
 }

void BTreeIndex::updateRootNode(PageKeyPair<int> keyPage, PageId currPage) {
	Page* root;
	PageId rootId; // new root page number
	bufMgr -> allocPage(file, rootId, root);
	// convert to proper structure
	NonLeafNodeInt *rootPage = (NonLeafNodeInt*) root;

	// Need update level ?
	rootPage -> keyArray[0] = keyPage.key;
	rootPage -> pageNoArray[0] = currPage;
	rootPage -> pageNoArray[1] = keyPage.pageNo;

	Page *headerInfo;
	bufMgr -> readPage(file, headerPageNum, headerInfo); // Reads the first page (meta)
	IndexMetaInfo *headerPage = (IndexMetaInfo *) headerInfo;
	rootPageNum = rootId;
	headerPage -> rootPageNo = rootId;

	// Update pins needed?
	bufMgr -> unPinPage(file, rootId, true);
}


// -----------------------------------------------------------------------------
// BTreeIndex::startScan
// -----------------------------------------------------------------------------

void BTreeIndex::startScan(const void* lowValParm,
				   const Operator lowOpParm,
				   const void* highValParm,
				   const Operator highOpParm)
{
    // Add your code below. Please do not remove this line.

	// end other scans
	if (scanExecuting) {
		endScan();
	}

	// initial conditions and exceptions
	if (lowOpParm != GT && lowOpParm != GTE) {
		throw BadOpcodesException();
	}
	if (highOpParm != LT && highOpParm != LTE) {
		throw BadOpcodesException();
	}
	
	lowValInt = *((int *)lowValParm);
	highValInt = *((int *)highValParm);
	if (lowValInt > highValInt) {
		throw BadScanrangeException();
	}

	lowOp = lowOpParm;
	highOp = highOpParm;

	// start new scan
	scanExecuting = true;
	nextEntry = startScanHelper(rootPageNum);

}

const int BTreeIndex::startScanHelper(PageId pageNum)
{
	NonLeafNodeInt* currNode;
	LeafNodeInt* child;

	currentPageNum = pageNum;
	bufMgr -> readPage(file, currentPageNum, currentPageData);
	currNode = (NonLeafNodeInt*)currentPageData;
	for (int i = 0; i < INTARRAYNONLEAFSIZE + 1; ++i) {
		if (i == INTARRAYNONLEAFSIZE || currNode -> pageNoArray[i + 1] == INVALID_NUMBER || currNode -> keyArray[i] > lowValInt) {
            // node is directly above leaf node 
            if (currNode -> level == 1) { 
                // read child page (leaf), then update
                currentPageNum = currNode -> pageNoArray[i];
                bufMgr -> readPage(file, currentPageNum, currentPageData);
                child = (LeafNodeInt*)currentPageData;
                // scan page, return index
                for(int j = 0; j < INTARRAYLEAFSIZE; ++j){
                    if (child -> ridArray[j] == INVALID_SLOT) {
                        break;
                    }
                    if ((lowOp == GT && child -> keyArray[j] > lowValInt)
                    || (lowOp == GTE && child -> keyArray[j] >= lowValInt)){
                        bufMgr -> unPinPage(file, pageNum, false);
                        return j;
                    }
                }   
                // if not found by the end
                bufMgr -> unPinPage(file, pageNum, false); 
                endScan();  
                throw NoSuchKeyFoundException();      
            } else { 
				// recursive call for nonLeaf nodes
                try {
                    int next = startScanHelper(currNode -> pageNoArray[i]);
                    bufMgr -> unPinPage(file, pageNum, false);
                    return next;
                } catch (NoSuchKeyFoundException) {
                    bufMgr -> unPinPage(file, pageNum, false);
                    throw NoSuchKeyFoundException(); 
                }
            }
        }
    }

	// impossible
    assert(false);
    return -1;

}


// -----------------------------------------------------------------------------
// BTreeIndex::scanNext
// -----------------------------------------------------------------------------

void BTreeIndex::scanNext(RecordId& outRid) 
{
    // Add your code below. Please do not remove this line.
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//
void BTreeIndex::endScan() 
{
    // Add your code below. Please do not remove this line.
}

}
