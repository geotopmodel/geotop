#define BOOST_TEST_DYN_LINK
#ifdef STAND_ALONE
#   define BOOST_TEST_MODULE Main
#endif
#include <boost/test/unit_test.hpp>

#include <AHQueue.h>
#include <boost/thread.hpp>
#include <boost/cstdint.hpp>

typedef adhoc::tool1::AHQueue<boost::uint32_t> uint32_q_t;

struct ToolsFixture
{

    uint32_q_t uint32Queue ;

    ToolsFixture()
    {
        BOOST_TEST_MESSAGE("setup tools");
    }

    ~ToolsFixture()
    {
    }
};

void queue_push_worker(const size_t pSize, uint32_q_t *pQ ) {
    for ( size_t i = 0 ; i < pSize ; i++ ) {
        pQ->push (i);
    }
    BOOST_TEST_MESSAGE("Push(" << pSize << ") Size(" << pQ->getSize() << ")");
}

void queue_wait_and_pop_worker(const size_t pSize, uint32_q_t *pQ ) {
    boost::uint32_t value ;
    size_t lPop = 0 ;
    while( not pQ->isEmpty () || lPop <= pSize )
    {
        if ( lPop == pSize )
            break ;

        pQ->waitAndPop(value);
        lPop += 1;
    }
    BOOST_TEST_MESSAGE("(wait) Pop(" << lPop << ") Size(" << pQ->getSize() << ")");
}

void queue_timed_wait_and_pop_worker(const size_t pSize, uint32_q_t *pQ ) {
    boost::uint32_t value ;
    size_t lPop = 0 ;
    while ( not pQ->isEmpty () || lPop <= pSize ) {
        if ( lPop == pSize )
            break ;

        if( pQ->waitAndPop(value, 5) )
            lPop += 1;
    }
    BOOST_TEST_MESSAGE("(timed_wait) Pop(" << lPop << ") Size(" << pQ->getSize() << ")");
}

void queue_pop_worker(const size_t pSize, uint32_q_t *pQ ) {
    boost::uint32_t value ;
    size_t lPop = 0 ;
    while ( not pQ->isEmpty () || lPop <= pSize ) {
        if ( lPop == pSize )
            break ;

        if( pQ->tryPop (value) )
            lPop += 1 ;
    }

    BOOST_TEST_MESSAGE("(try) Pop(" << lPop << ") Size(" << pQ->getSize() << ")");
}

BOOST_FIXTURE_TEST_SUITE(test_queue, ToolsFixture)

BOOST_AUTO_TEST_CASE(enqueue_and_dequeue)
{
    // concurrent push test
    {
        size_t s0 = 100000 ;
        size_t s1 = 100000 ;
        size_t s2 = 100000 ;
        size_t s3 = 100000 ;

        boost::thread w0 (queue_push_worker, s0, &uint32Queue) ;
        boost::thread w1 (queue_push_worker, s1, &uint32Queue) ;
        boost::thread w2 (queue_push_worker, s2, &uint32Queue) ;
        boost::thread w3 (queue_push_worker, s3, &uint32Queue) ;

        w0.join ();
        w1.join ();
        w2.join ();
        w3.join ();

        size_t num_of_elements = s0 + s1 + s2 + s3 ;

        BOOST_CHECK( uint32Queue.getSize () == num_of_elements );
    }

    // concurrent try pop test
    {
        size_t s0 = 100000 ;
        size_t s1 = 100000 ;
        size_t s2 = 100000 ;
        size_t s3 = 100000 ;

        boost::thread w0 (queue_pop_worker, s0, &uint32Queue) ;
        boost::thread w1 (queue_pop_worker, s1, &uint32Queue) ;
        boost::thread w2 (queue_pop_worker, s2, &uint32Queue) ;
        boost::thread w3 (queue_pop_worker, s3, &uint32Queue) ;

        w0.join ();
        w1.join ();
        w2.join ();
        w3.join ();

        BOOST_CHECK( uint32Queue.getSize () == 0 );
    }
}

BOOST_AUTO_TEST_CASE(concurrent_push_trypop)
{
    //concurrent push/trypop
    {
        size_t s0 = 100000 ;
        size_t s1 = 100000 ;
        size_t s2 = 100000 ;
        size_t s3 = 100000 ;

        boost::thread w0 (queue_push_worker, s0, &uint32Queue) ;
        boost::thread w1 (queue_push_worker, s1, &uint32Queue) ;
        boost::thread w2 (queue_push_worker, s2, &uint32Queue) ;
        boost::thread w3 (queue_push_worker, s3, &uint32Queue) ;
        boost::thread w4 (queue_pop_worker, s0, &uint32Queue) ;
        boost::thread w5 (queue_pop_worker, s1, &uint32Queue) ;
        boost::thread w6 (queue_pop_worker, s2, &uint32Queue) ;
        boost::thread w7 (queue_pop_worker, s3, &uint32Queue) ;

        w0.join ();
        w1.join ();
        w2.join ();
        w3.join ();
        w4.join ();
        w5.join ();
        w6.join ();
        w7.join ();

        BOOST_CHECK( uint32Queue.getSize () == 0 );
    }
}

BOOST_AUTO_TEST_CASE(concurrent_push_and_waitpop_test)
{
    //concurrent push and waitpop test
    {
        size_t s0 = 100000 ;
        size_t s1 = 100000 ;
        size_t s2 = 100000 ;
        size_t s3 = 100000 ;

        boost::thread w0 (queue_push_worker, s0, &uint32Queue) ;
        boost::thread w1 (queue_push_worker, s1, &uint32Queue) ;
        boost::thread w2 (queue_push_worker, s2, &uint32Queue) ;
        boost::thread w3 (queue_push_worker, s3, &uint32Queue) ;
        boost::thread w4 (queue_wait_and_pop_worker, s0, &uint32Queue) ;
        boost::thread w5 (queue_wait_and_pop_worker, s1, &uint32Queue) ;
        boost::thread w6 (queue_wait_and_pop_worker, s2, &uint32Queue) ;
        boost::thread w7 (queue_wait_and_pop_worker, s3, &uint32Queue) ;

        w0.join ();
        w1.join ();
        w2.join ();
        w3.join ();
        w4.join ();
        w5.join ();
        w6.join ();
        w7.join ();

        BOOST_CHECK( uint32Queue.getSize () == 0 );
    }
}

BOOST_AUTO_TEST_CASE(concurrent_push_and_timed_waitpop_test)
{
    //concurrent push and timed waitpop test
    {
        size_t s0 = 100000 ;
        size_t s1 = 100000 ;
        size_t s2 = 100000 ;
        size_t s3 = 100000 ;

        boost::thread w0 (queue_push_worker, s0, &uint32Queue) ;
        boost::thread w1 (queue_push_worker, s1, &uint32Queue) ;
        boost::thread w2 (queue_push_worker, s2, &uint32Queue) ;
        boost::thread w3 (queue_push_worker, s3, &uint32Queue) ;
        boost::thread w4 (queue_timed_wait_and_pop_worker, s0, &uint32Queue) ;
        boost::thread w5 (queue_timed_wait_and_pop_worker, s1, &uint32Queue) ;
        boost::thread w6 (queue_timed_wait_and_pop_worker, s2, &uint32Queue) ;
        boost::thread w7 (queue_timed_wait_and_pop_worker, s3, &uint32Queue) ;

        w0.join ();
        w1.join ();
        w2.join ();
        w3.join ();
        w4.join ();
        w5.join ();
        w6.join ();
        w7.join ();

        BOOST_CHECK( uint32Queue.getSize () == 0 );
    }
}

BOOST_AUTO_TEST_SUITE_END()

