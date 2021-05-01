/**
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
	this->attrByteOffset = attrByteOffset;
	attributeType = attrType;
	scanExecuting = false;
	leafOccupancy = INTARRAYLEAFSIZE;
	nodeOccupancy = INTARRAYNONLEAFSIZE;

	std::ostringstream idxStr;
	idxStr << relationName << "." << attrByteOffset;
	outIndexName = idxStr.str();

	Page *headerPage, *rootPage;
	IndexMetaInfo* meta;

	try {
		file = new BlobFile(outIndexName, false);

		headerPageNum = file -> getFirstPageNo();
		bufMgr -> readPage(file, headerPageNum, headerPage);
		meta = (IndexMetaInfo*)headerPage;
		rootPageNum = meta -> rootPageNo; // Held in the meta page

		bufMgr -> unPinPage(file, headerPageNum, false); 
	}
	catch (FileNotFoundException e) {
		file = new BlobFile(outIndexName, true);
		// Allocate header page
		bufMgr -> allocPage(file, headerPageNum, headerPage);
		memset(headerPage, 0, Page::SIZE);
		// Fills in meta info
		meta = (IndexMetaInfo*) headerPage;
		// Allocate root page
    	bufMgr -> allocPage(file, rootPageNum, (Page *&)rootPage); 
		memset(rootPage, 0, Page::SIZE);

		strcpy(meta -> relationName, relationName.c_str());
		meta -> attrByteOffset = attrByteOffset;
		meta -> attrType = attrType;
		meta -> rootPageNo = rootPageNum;
		initialRootPageNum = rootPageNum;
		
		LeafNodeInt *root = (LeafNodeInt *)rootPage;
    	root -> rightSibPageNo = 0;

		bufMgr -> unPinPage(file, headerPageNum, true); 
		bufMgr -> unPinPage(file, rootPageNum, true); 

		// Scans Records
		FileScan scan = FileScan(relationName, bufMgr);
		
		try {
			RecordId rid;
			while(1) {
				scan.scanNext(rid);
				std::string fileRecord = scan.getRecord();
				const char *record = fileRecord.c_str();
				insertEntry(record + attrByteOffset, rid);
			}
		}
		catch (EndOfFileException e) {
			
		}
	}
}


// -----------------------------------------------------------------------------
// BTreeIndex::~BTreeIndex -- destructor
// -----------------------------------------------------------------------------

BTreeIndex::~BTreeIndex()
{
	/// stop the ongoing scanning
	scanExecuting = false;
	bufMgr -> flushFile(BTreeIndex::file);
	delete file;
	file = nullptr;
}

/**
* Insert a key and a record ID to the appropriate position and array of the given leaf node.
*
* @param key   key to be inserted to keyArray
* @param rid   rid to be inserted to ridArray
* @param node  corresponding leaf node to be modified
* @param index index of a particular key and rid to be inserted
*/
void BTreeIndex::insertLeafNode(int key, RecordId rid, LeafNodeInt* node, int index) {
	int i = INTARRAYLEAFSIZE; // Number of keys in the leaf node // Maybe change to leafOccupancy ?
	while((node -> ridArray[i-1].page_number == 0) && i > 0) { // Gets to the end of leafNode
		i--;
	}
	while((node -> keyArray[i-1] > key) && i > 0) { // Shifts the previous ridKey pairs
		node -> keyArray[i] = node -> keyArray[i-1];
		node -> ridArray[i] = node -> ridArray[i-1];
		i--;
	}
	node -> keyArray[i] = key; // Inserts ridKey to the leafNode
	node -> ridArray[i] = rid;

	// Increment size
	node->size += 1;
}

/**
* A helper method for inserting a nonleaf node in the B Tree
*
* @param key   inserted key
* @param childPageId   inserted PageId
* @param node  node that is modified
* @param index index of a particular key and pageid to be inserted
*/
void BTreeIndex::insertNonLeafNode(int key, PageId childPageId, NonLeafNodeInt* node, int index) {
	int i = nodeOccupancy;
	while((node -> pageNoArray[i] == 0) && i >= 0) {
		i--;
	}
	while((node -> keyArray[i-1] > key) && i > 0) {
		node -> keyArray[i] = node -> keyArray[i-1]; // Shifts the previous keyPage pairs
		node -> pageNoArray[i+1] = node -> pageNoArray[i];
		i--;
	}

	node -> keyArray[i] = key; // Inserts keyPage pair to the node
	node -> pageNoArray[i+1] = childPageId;

	node -> size += 1; // Increase node size
}

/**
* Helper function to do the recurion of inserting a key and a rid to the tree.
*
* @param key   key to be inserted
* @param rid   rid to be inserted
* @param currPageId  the current page which holds the node
* @param midVal will pass back the middle value to its parent
* @param addedPage will pass back the splitted node's pageID to its parent
* @param isLeaf whether node is a leaf
*/
void BTreeIndex::insertEntryHelper(const int key,
					const RecordId rid,
					PageId currPageId,
					int* midVal,
					PageId* addedPage,
					bool isLeaf) 
{
	Page *currNode;
  	bufMgr->readPage(file, currPageId, currNode);

  	// urrent node is leaf 
  	if (isLeaf) {
  		LeafNodeInt* currLeafNode = (LeafNodeInt *)currNode;

  		// index of new key
  		int index = currLeafNode->size;
  		for (int i = 0; i < currLeafNode->size; i++) {
  			if (currLeafNode->keyArray[i] > key) {
  				index = i;
  				break;
  			}
  		}

  		// insert to the leaf
  		if (currLeafNode->ridArray[INTARRAYLEAFSIZE - 1].page_number == Page::INVALID_NUMBER) {
			insertLeafNode(key, rid, currLeafNode, index);
  			bufMgr->unPinPage(file, currPageId, true);
  			*midVal = 0;
  			*addedPage = 0;
  		} 
  		// else, split
  		else {
  			splitLeaf(currLeafNode, currPageId, key, rid, index, midVal, addedPage);
		}
  	}
  	// current node is not a leaf node
  	else {
  		NonLeafNodeInt *currNonLeaf = (NonLeafNodeInt *)currNode;

  		// ind the  child 
  		int childIndex = currNonLeaf->size;
  		for (int i = 0; i < currNonLeaf->size; i++) {
  			if (currNonLeaf->keyArray[i] > key) {
  				childIndex = i;
  				break;
  			}
  		}
  		PageId currChildId = currNonLeaf->pageNoArray[childIndex];

  		// recurse to child
		int childMidKey;
		PageId newChildId;
		insertEntryHelper(key, rid, currChildId, &childMidKey, &newChildId, currNonLeaf->level == 1);

		// no split
		if ((int) newChildId == 0) {
		  bufMgr->unPinPage(file, currPageId, false);
		  *midVal = 0;
		  *addedPage = 0;
		}
		// child node exists
		else {
	  		int index = currNonLeaf->size;
	  		for (int i = 0; i < currNonLeaf->size; i++) {
	  			if (currNonLeaf->keyArray[i] > childMidKey) {
	  				index = i;
	  				break;
	  			}
	  		}

	  		// nsert middle key to the node
  			if (currNonLeaf->pageNoArray[INTARRAYNONLEAFSIZE] == Page::INVALID_NUMBER) {
  				insertNonLeafNode(childMidKey, newChildId, currNonLeaf, index);
  				bufMgr->unPinPage(file, currPageId, true);
	  			*midVal = 0;
	  			*addedPage = 0;
  			}
  			// else split
  			else {
	  			PageId newPageId;
	  			Page* newNode;
	  			bufMgr->allocPage(file, newPageId, newNode);
	  			memset(newNode, 0, Page::SIZE);
	  			NonLeafNodeInt* newNonLeaf = (NonLeafNodeInt *)newNode;
	  			newNonLeaf->level = currNonLeaf->level;

	  			// splitting
	  			int midIndex = INTARRAYNONLEAFSIZE / 2;
	  			// new child will be pushed
	  			if (index == midIndex) {
	  				for (int i = midIndex; i < INTARRAYNONLEAFSIZE; i++) {
						newNonLeaf->keyArray[i-midIndex] = currNonLeaf->keyArray[i];
						newNonLeaf->pageNoArray[i-midIndex+1] = currNonLeaf->pageNoArray[i+1];
						currNonLeaf->keyArray[i] = 0;
						currNonLeaf->pageNoArray[i+1] = Page::INVALID_NUMBER;
					}
					newNonLeaf->pageNoArray[0] = newChildId;
					currNonLeaf->size = midIndex;
					newNonLeaf->size = INTARRAYNONLEAFSIZE - midIndex;

		  			bufMgr->unPinPage(file, currPageId, true);
		  			bufMgr->unPinPage(file, newPageId, true);

					// return new values 
		  			*midVal = childMidKey;
		  			*addedPage = newPageId;
	  			}

	  			else {
	  				if (INTARRAYNONLEAFSIZE % 2 == 0 && index < midIndex) {
	  					midIndex -= 1;
	  				}
	  				for (int i = midIndex + 1; i < INTARRAYNONLEAFSIZE; i++) {
						newNonLeaf->keyArray[i-midIndex-1] = currNonLeaf->keyArray[i];
						newNonLeaf->pageNoArray[i-midIndex-1] = currNonLeaf->pageNoArray[i];
						currNonLeaf->keyArray[i] = 0;
						currNonLeaf->pageNoArray[i] = Page::INVALID_NUMBER;
					}

					// return new values, clear, set new size
					*midVal = currNonLeaf->keyArray[midIndex];
					*addedPage = newPageId;
					currNonLeaf->keyArray[midIndex] = 0;
					currNonLeaf->size = midIndex;
					newNonLeaf->size = INTARRAYNONLEAFSIZE - midIndex - 1;

					// insert value to left or right node 
					if (index < INTARRAYNONLEAFSIZE / 2) {
						insertNonLeafNode(childMidKey, newChildId, currNonLeaf, index);
					} else {
						insertNonLeafNode(childMidKey, newChildId, newNonLeaf, index - midIndex);
					}

		  			bufMgr->unPinPage(file, currPageId, true);
		  			bufMgr->unPinPage(file, newPageId, true);
	  			}
  			}
		}
  	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::insertEntry
// -----------------------------------------------------------------------------

void BTreeIndex::insertEntry(const void *key, const RecordId rid) 
{
	PageId newPageId;
	int midValChild;	
	insertEntryHelper(*(int*)key, rid, rootPageNum, &midValChild, &newPageId, rootIsLeaf); // Gets new midVal and PageId in a split

	if ((int) newPageId != 0) { // Checks for split in root
		Page* root; // New root page number
	  	PageId rootId;
		bufMgr->allocPage(file, rootId, root);
		memset(root, 0, Page::SIZE);
		NonLeafNodeInt* rootPage = (NonLeafNodeInt *)root; // Converts to proper structure
		
		rootPage -> size = 1; // Updates the root page
		rootPage -> keyArray[0] = midValChild;
		rootPage -> pageNoArray[0] = rootPageNum;
		rootPage -> pageNoArray[1] = newPageId;

		if (!rootIsLeaf) {
			rootPage -> level = 0;
		} 
		else {
			rootPage -> level = 1;
		}
		rootIsLeaf = false;

		Page *newMeta; // Updates the first page
		bufMgr -> readPage(file, headerPageNum, newMeta); // Reads the first page (meta)
		IndexMetaInfo *metaPage = (IndexMetaInfo*) newMeta;
		metaPage -> rootPageNo = newPageId;
		rootPageNum = newPageId;

		bufMgr -> unPinPage(file, headerPageNum, true); // Unpins at the end of inserting
		bufMgr -> unPinPage(file, newPageId, true);
	}
}

const void BTreeIndex::splitLeaf(LeafNodeInt *currLeafNode, PageId leafPageNum, const int key, const RecordId rid, const int index, int* midVal, PageId* addedPage){
	PageId newPageId;
	Page* newNode;
	bufMgr->allocPage(file, newPageId, newNode);
	memset(newNode, 0, Page::SIZE);
	LeafNodeInt* newLeafNode = (LeafNodeInt *)newNode;

	int midIndex = INTARRAYLEAFSIZE / 2;
	if (INTARRAYLEAFSIZE % 2 == 1 && index > midIndex) {
		midIndex = midIndex + 1;
	}

	for (int i = midIndex; i < INTARRAYLEAFSIZE; i++) {
		newLeafNode->keyArray[i-midIndex] = currLeafNode->keyArray[i];
		newLeafNode->ridArray[i-midIndex] = currLeafNode->ridArray[i];
		currLeafNode->keyArray[i] = 0;
		currLeafNode->ridArray[i].page_number = Page::INVALID_NUMBER;
	}

	newLeafNode->size = INTARRAYLEAFSIZE - midIndex;
	currLeafNode->size = midIndex;

	// insert to right 
	if (index > INTARRAYLEAFSIZE / 2) {
		insertLeafNode(key, rid, newLeafNode, index - midIndex);
	}
	// insert to left 
	else {
		insertLeafNode(key, rid, currLeafNode, index);
	}

	newLeafNode->rightSibPageNo = currLeafNode->rightSibPageNo;
	currLeafNode->rightSibPageNo = newPageId;

	bufMgr->unPinPage(file, leafPageNum, true);
	bufMgr->unPinPage(file, newPageId, true);

	*midVal = newLeafNode->keyArray[0];
	*addedPage = newPageId;	
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

/**
 * Helper method for the startScan method to update nextEntry number
 * @param pageNum 	 
 * @return 	int 
 */
const int BTreeIndex::startScanHelper(PageId pageNum)
{
	NonLeafNodeInt* currNode;
	LeafNodeInt* child;

	currentPageNum = pageNum;
	bufMgr -> readPage(file, currentPageNum, currentPageData);
	currNode = (NonLeafNodeInt*)currentPageData;
	for (int i = 0; i < INTARRAYNONLEAFSIZE + 1; ++i) {
		if (i == INTARRAYNONLEAFSIZE || currNode -> pageNoArray[i + 1] == Page::INVALID_NUMBER || currNode -> keyArray[i] > lowValInt) {
            // node is directly above leaf node 
            if (currNode -> level == 1) { 
                // read child page (leaf), then update
                currentPageNum = currNode -> pageNoArray[i];
                bufMgr -> readPage(file, currentPageNum, currentPageData);
                child = (LeafNodeInt*)currentPageData;
                // scan page, return index
                for (int j = 0; j < INTARRAYLEAFSIZE; ++j) {
                    if (child -> ridArray[j].page_number == Page::INVALID_NUMBER) {
                        break;
                    }
                    if ((lowOp == GT && child -> keyArray[j] > lowValInt)
                    || (lowOp == GTE && child -> keyArray[j] >= lowValInt)) {
                        bufMgr -> unPinPage(file, pageNum, false);
                        return j;
                    }
                }   
                // if not found by the end
                bufMgr -> unPinPage(file, pageNum, false); 
                endScan();  
                throw NoSuchKeyFoundException();      
            } 
			else { 
				// recursive call for nonLeaf nodes
                try {
                    int next = startScanHelper(currNode -> pageNoArray[i]);
                    bufMgr -> unPinPage(file, pageNum, false);
                    return next;
                } catch (NoSuchKeyFoundException e) {
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
/**
 * Fetch the record id of the next index entry that matches the scan.
 * Return the next record from current page being scanned. If current page has been scanned to its entirety, 
 * move on to the right sibling of current page, if any exists, to start scanning that page. Make sure to unpin any pages that are no longer required.
 * @param outRid	RecordId of next record found that satisfies the scan criteria returned in this
 * @throws ScanNotInitializedException If no scan has been initialized.
 * @throws IndexScanCompletedException If no more records, satisfying the scan criteria, are left to be scanned.
 **/
void BTreeIndex::scanNext(RecordId& outRid) 
{
    // Add your code below. Please do not remove this line.
	if (scanExecuting == false) {
		throw ScanNotInitializedException();
	}

	LeafNodeInt* leafNode = (LeafNodeInt*) currentPageData;
	int nextIndex;
	if (leafNode->ridArray[nextEntry].page_number == Page::INVALID_NUMBER) {
		throw IndexScanCompletedException();
	}
	// done scanning
	else if (nextEntry == INT_MAX) {
		throw IndexScanCompletedException();
	}
	else {
		nextIndex = nextEntry;
	}

	if ((leafNode -> keyArray[nextIndex] <= highValInt) && highOp == LTE) {
		outRid = leafNode -> ridArray[nextIndex];
	}
	else if ((leafNode -> keyArray[nextIndex] < highValInt) && highOp == LT) {
		outRid = leafNode -> ridArray[nextIndex];
	}
	// done scanning
	else {
		throw IndexScanCompletedException(); 
	}

	if (nextEntry < INTARRAYLEAFSIZE - 1) {
		if (leafNode -> ridArray[nextEntry + 1].page_number != Page::INVALID_NUMBER) {
			nextEntry++;
		}
		else {
			if (leafNode -> rightSibPageNo == Page::INVALID_NUMBER) {
				nextEntry = INT_MAX;
			}
			else { 
				PageId nextPageId = currentPageNum;
				currentPageNum = leafNode -> rightSibPageNo; 
				nextEntry = 0;
				bufMgr -> readPage(file, currentPageNum, currentPageData); 
				bufMgr -> unPinPage(file, nextPageId, false); 
			}
		}
	}
}

// -----------------------------------------------------------------------------
// BTreeIndex::endScan
// -----------------------------------------------------------------------------
//
/**
 * Terminate the current scan. Unpin any pinned pages. Reset scan specific variables.
 * @throws ScanNotInitializedException If no scan has been initialized.
**/
void BTreeIndex::endScan() 
{
    // Add your code below. Please do not remove this line.
	//check if scan has been initialized
	if (scanExecuting == false){
		throw ScanNotInitializedException();
	}
	//terminate scan
	scanExecuting = false;
	//unpin pages
	bufMgr -> unPinPage(file, currentPageNum, false);
	//reset scan specific variables
	nextEntry = -1;
	currentPageNum = -1;
	currentPageData = nullptr;
	lowValInt = -1;
	highValInt = -1;
	lowOp = (Operator)-1;
	highOp = (Operator)-1;
}

}
